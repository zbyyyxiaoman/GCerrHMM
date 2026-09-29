#!/usr/bin/env python3
"""Apply the pre-registered GC-bin ablation rule (docs/gc_claim_evidence_chain.md).

Rule: "1-bin is worst" is accepted only if, **in both species**, the 1-bin
model has the lowest mean GC-fidelity Pearson r and the highest mean MAD, and
the gap exceeds the pooled standard deviation over the three seeds. Anything
else selects the contracted claim, and the contracted branch must then be
written with the same prominence as the accepted one.

usage: python3 ablation_claim_decision.py --project-dir ~/errhmm_project
"""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from pathlib import Path

SPECIES = [("Ecoli", "E. coli"), ("Athaliana", "A. thaliana")]
BINS = [1, 5, 10, 20]
SEEDS = [11, 22, 33]


def load(path: Path) -> dict:
    if not path.exists() or path.stat().st_size == 0:
        return {}
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError:
        return {}


def mean_sd(values: list[float]) -> tuple[float | None, float | None]:
    values = [v for v in values if v is not None]
    if not values:
        return None, None
    if len(values) == 1:
        return values[0], 0.0
    return statistics.mean(values), statistics.stdev(values)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", required=True)
    args = parser.parse_args()

    project = Path(args.project_dir).expanduser()
    fid_dir = project / "results" / "gc_fidelity"
    out_dir = project / "results" / "stats"
    out_dir.mkdir(parents=True, exist_ok=True)

    table = {}
    for key, label in SPECIES:
        for bins in BINS:
            r_values, mad_values, seeds_used = [], [], []
            for seed in SEEDS:
                path = fid_dir / f"ablation_{key}_{bins}bin_seed{seed}.gc.json"
                payload = load(path)
                if not payload:
                    continue
                if payload.get("pearson_r") is not None:
                    r_values.append(float(payload["pearson_r"]))
                if payload.get("mad") is not None:
                    mad_values.append(float(payload["mad"]))
                seeds_used.append(seed)
            r_mean, r_sd = mean_sd(r_values)
            mad_mean, mad_sd = mean_sd(mad_values)
            table[(key, bins)] = {
                "species": label,
                "gc_bins": bins,
                "n": len(seeds_used),
                "seeds": seeds_used,
                "r_mean": r_mean,
                "r_sd": r_sd,
                "mad_mean": mad_mean,
                "mad_sd": mad_sd,
            }

    complete = all(table[(key, bins)]["n"] == 3
                   for key, _ in SPECIES for bins in BINS)

    rows = list(table.values())
    csv_path = out_dir / "gc_bins_claim_decision.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    # ---- apply the rule -------------------------------------------------
    verdicts = {}
    for key, label in SPECIES:
        entries = {b: table[(key, b)] for b in BINS}
        r1, r1sd = entries[1]["r_mean"], entries[1]["r_sd"]
        mad1, mad1sd = entries[1]["mad_mean"], entries[1]["mad_sd"]
        if r1 is None or mad1 is None:
            verdicts[key] = {"worst": None, "reason": "no 1-bin measurement"}
            continue
        rivals = [b for b in BINS if b != 1 and entries[b]["r_mean"] is not None]
        r_worst = all(r1 < entries[b]["r_mean"] for b in rivals) if rivals else None
        mad_worst = all(mad1 > entries[b]["mad_mean"] for b in rivals) if rivals else None
        # gap must exceed the pooled SD of the pair being compared
        r_significant, mad_significant = True, True
        for b in rivals:
            r_gap = entries[b]["r_mean"] - r1
            pooled_r = ((r1sd or 0) ** 2 + (entries[b]["r_sd"] or 0) ** 2) ** 0.5
            if r_gap <= pooled_r:
                r_significant = False
            mad_gap = mad1 - entries[b]["mad_mean"]
            pooled_m = ((mad1sd or 0) ** 2 + (entries[b]["mad_sd"] or 0) ** 2) ** 0.5
            if mad_gap <= pooled_m:
                mad_significant = False
        worst = bool(r_worst and mad_worst and r_significant and mad_significant)
        verdicts[key] = {
            "worst": worst,
            "r_worst": r_worst,
            "mad_worst": mad_worst,
            "r_gap_exceeds_sd": r_significant,
            "mad_gap_exceeds_sd": mad_significant,
        }

    accepted = complete and all(v.get("worst") for v in verdicts.values())

    lines = [
        "# GC-bin ablation: pre-registered claim decision",
        "",
        "Rule: `docs/gc_claim_evidence_chain.md`. Primary instrument: "
        "GC-conditional error fidelity (Pearson r against the real per-bin "
        "error rate, and MAD).",
        "",
        f"Data completeness: {'complete' if complete else 'INCOMPLETE - verdict provisional'}",
        "",
        "| species | GC bins | n | mean r | r SD | mean MAD (pp) | MAD SD (pp) |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        r = "-" if row["r_mean"] is None else f"{row['r_mean']:.3f}"
        rsd = "-" if row["r_sd"] is None else f"{row['r_sd']:.3f}"
        mad = "-" if row["mad_mean"] is None else f"{row['mad_mean'] * 100:.3f}"
        madsd = "-" if row["mad_sd"] is None else f"{row['mad_sd'] * 100:.3f}"
        lines.append(f"| {row['species']} | {row['gc_bins']} | {row['n']} | "
                     f"{r} | {rsd} | {mad} | {madsd} |")

    lines += ["", "## Per-species checks", ""]
    for key, label in SPECIES:
        v = verdicts[key]
        lines.append(
            f"* **{label}**: 1-bin worst = {v.get('worst')} "
            f"(lowest r: {v.get('r_worst')}, highest MAD: {v.get('mad_worst')}, "
            f"r gap > SD: {v.get('r_gap_exceeds_sd')}, "
            f"MAD gap > SD: {v.get('mad_gap_exceeds_sd')})"
        )

    lines += ["", "## Branch", ""]
    if not complete:
        lines.append(
            "**PROVISIONAL**: the ablation is not complete yet, so no branch "
            "is selected. Re-run this script once every species/bin/seed point "
            "has a GC-fidelity result; the decision must not be read out of a "
            "partial table."
        )
    elif accepted:
        lines.append(
            "**Accepted branch**: GC-conditioned error modelling improves "
            "GC-conditional error fidelity. The composite panel is demoted to "
            "supporting evidence and the claim is stated in terms of the "
            "fidelity instrument."
        )
    else:
        lines.append(
            "**Contracted branch**: the claim contracts to a trainable "
            "GC-aware simulation framework with competitive compositional and "
            "downstream fidelity; GC conditioning is reported as a design "
            "feature rather than a performance claim. Write this branch with "
            "the same prominence as the accepted one would have received."
        )

    md_path = out_dir / "gc_bins_claim_decision.md"
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print(f"\nwritten: {csv_path}")
    print(f"written: {md_path}")


if __name__ == "__main__":
    main()
