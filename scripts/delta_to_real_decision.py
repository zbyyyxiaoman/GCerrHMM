#!/usr/bin/env python3
"""Apply the pre-registered delta-to-real rule to the chr21 30x panel.

The rule lives in docs/delta_to_real_framing.md and was written down before the
coverage-matched real anchor existed:

  R3 alignment : |base_identity(tool) - base_identity(real)|
  R4 variant   : |f1_score(tool)      - f1_score(real)|
  R5 assembly  : Borda rank of |delta identity| and |log2(N50 ratio)|

GCerrHMM needs a top-2 finish in at least two of the three layers for the
"downstream consistency" headline; otherwise the paper falls back to the
compositional-fidelity headline and credits PBSim3-sample's sampling fidelity.

usage: python3 delta_to_real_decision.py --project-dir ~/errhmm_project
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

from common_io import read_json, sha256_file

SPECIES = "Hsapiens_chr21"
SUBSPECIES = "30x"
TOOLS = [
    ("errhmm", "GCerrHMM"),
    ("nanosim", "NanoSim"),
    ("pbsim", "PBSim3-sample"),
    ("badread", "badread"),
]
REAL = "real_ont30"


def table(dirs, name: str) -> dict:
    for directory in dirs:
        payload = read_json(directory / name, default={}, ignore_errors=True)
        if payload:
            return payload
    return {}


def rank(values: dict[str, float]) -> dict[str, float]:
    """Average ranks for ties, so no tool gets a favourable input order."""
    ordered = sorted(values.items(), key=lambda kv: kv[1])
    ranks: dict[str, float] = {}
    index = 0
    while index < len(ordered):
        end = index + 1
        while end < len(ordered) and ordered[end][1] == ordered[index][1]:
            end += 1
        average_rank = (index + 1 + end) / 2
        for key, _ in ordered[index:end]:
            ranks[key] = average_rank
        index = end
    return ranks


def fmt(value, spec=".4f"):
    return format(value, spec) if isinstance(value, (int, float)) else "-"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", required=True)
    parser.add_argument(
        "--exclude-layer", action="append", default=[],
        choices=["R3", "R4", "R5"],
        help="Layer whose real anchor is structurally unavailable (for example "
             "R4 when the truth set is a simulation truth); it is printed but "
             "not counted towards the >=2-of-3 rule.",
    )
    args = parser.parse_args()
    excluded = set(args.exclude_layer)

    project = Path(args.project_dir).expanduser()
    dirs = [project / "results" / "framework" / "tables",
            project / "results" / "tables"]
    out_dir = project / "results" / "framework" / "stats"
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for key, label in TOOLS:
        mapping = table(dirs, f"mapping_{key}_{SPECIES}_{SUBSPECIES}.json")
        variant = table(dirs, f"variant_{key}_{SPECIES}_{SUBSPECIES}.json")
        assembly = table(dirs, f"assembly_flye_{key}_{SPECIES}_{SUBSPECIES}.json")
        if not mapping or not variant or not assembly:
            raise SystemExit(
                f"30x panel table missing for {label} ({key})"
            )
        rows.append({
            "tool_key": key,
            "tool_label": label,
            "alignment_identity": mapping.get("base_identity"),
            "variant_f1": variant.get("f1_score"),
            "snp_f1": variant.get("snp_f1_score"),
            "indel_f1": variant.get("indel_f1_score"),
            "assembly_identity": assembly.get("reference_identity"),
            "assembly_n50": assembly.get("n50"),
            "assembly_contigs": assembly.get("num_contigs"),
        })

    real_mapping = table(dirs, f"mapping_{REAL}_{SPECIES}.json")
    real_variant = table(dirs, f"variant_{REAL}_{SPECIES}.json")
    real_assembly = table(dirs, f"assembly_flye_{REAL}_{SPECIES}.json")
    real = {
        "alignment_identity": real_mapping.get("base_identity"),
        "variant_f1": real_variant.get("f1_score"),
        "snp_f1": real_variant.get("snp_f1_score"),
        "indel_f1": real_variant.get("indel_f1_score"),
        "assembly_identity": real_assembly.get("reference_identity"),
        "assembly_n50": real_assembly.get("n50"),
        "assembly_contigs": real_assembly.get("num_contigs"),
    }

    missing = [name for name, value in real.items() if value is None]
    if missing:
        raise SystemExit(f"real anchor incomplete, missing: {missing}")

    layer_distances: dict[str, dict[str, float]] = {"R3": {}, "R4": {}, "R5": {}}
    ranks: dict[str, dict[str, float]] = {}

    for row in rows:
        row["delta_identity"] = abs(row["alignment_identity"] - real["alignment_identity"])
        row["delta_f1"] = abs(row["variant_f1"] - real["variant_f1"])
        row["delta_asm_identity"] = abs(row["assembly_identity"] - real["assembly_identity"])
        row["n50_ratio"] = (row["assembly_n50"] / real["assembly_n50"]
                            if row["assembly_n50"] and real["assembly_n50"] else None)
        layer_distances["R3"][row["tool_key"]] = row["delta_identity"]
        layer_distances["R4"][row["tool_key"]] = row["delta_f1"]

    # R5 is scale free: rank on each sub-metric, then average the two ranks.
    identity_rank = rank({r["tool_key"]: r["delta_asm_identity"] for r in rows})
    ratio_rank = rank({r["tool_key"]: abs(math.log2(r["n50_ratio"]))
                       if r["n50_ratio"] else float("inf") for r in rows})
    for row in rows:
        row["asm_rank_identity"] = identity_rank[row["tool_key"]]
        row["asm_rank_n50"] = ratio_rank[row["tool_key"]]
        row["asm_borda"] = (identity_rank[row["tool_key"]]
                            + ratio_rank[row["tool_key"]]) / 2
        layer_distances["R5"][row["tool_key"]] = row["asm_borda"]

    for layer in ("R3", "R4", "R5"):
        ranks[layer] = rank(layer_distances[layer])

    for row in rows:
        row["rank_R3"] = ranks["R3"][row["tool_key"]]
        row["rank_R4"] = ranks["R4"][row["tool_key"]]
        row["rank_R5"] = ranks["R5"][row["tool_key"]]
        row["top2_layers"] = sum(
            1 for layer in ("R3", "R4", "R5")
            if layer not in excluded and row[f"rank_{layer}"] <= 2
        )

    errhmm = next(r for r in rows if r["tool_key"] == "errhmm")
    counted = [layer for layer in ("R3", "R4", "R5") if layer not in excluded]
    headline_branch = errhmm["top2_layers"] >= min(2, len(counted))

    fields = list(rows[0].keys())
    fields += ["excluded_R3", "excluded_R4", "excluded_R5"]
    for row in rows:
        for layer in ("R3", "R4", "R5"):
            row[f"excluded_{layer}"] = layer in excluded
    csv_path = out_dir / "delta_to_real_decision.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "# chr21 30x delta-to-real decision",
        "",
        "Rule: `docs/delta_to_real_framing.md` (pre-registered 2026-09-25).",
        "Data: coverage-matched 30x panel "
        "(`results/framework/tables/*_Hsapiens_chr21_30x.json`).",
        "Tie policy: average ranks; top-2 requires an average rank <= 2.",
        "",
        "## Real anchor",
        "",
        f"- base identity: {fmt(real['alignment_identity'])}",
        f"- variant F1: {fmt(real['variant_f1'])} "
        f"(SNP {fmt(real['snp_f1'])}, indel {fmt(real['indel_f1'])})",
        f"- assembly: {real['assembly_contigs']} contigs, "
        f"N50 {fmt(real['assembly_n50'], '.0f')} bp, "
        f"identity {fmt(real['assembly_identity'])}",
        "",
        "## Per-tool distances",
        "",
        "| tool | R3 identity | ΔR3 | R4 F1 | ΔR4 | R5 identity | R5 N50 | "
        "rank R3 | rank R4 | rank R5 | top-2 layers |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        rank_display = {
            layer: (
                "excluded"
                if layer in excluded
                else fmt(row[f"rank_{layer}"], ".1f")
            )
            for layer in ("R3", "R4", "R5")
        }
        lines.append(
            f"| {row['tool_label']} | {fmt(row['alignment_identity'])} | "
            f"{fmt(row['delta_identity'])} | {fmt(row['variant_f1'])} | "
            f"{fmt(row['delta_f1'])} | {fmt(row['assembly_identity'])} | "
            f"{fmt(row['assembly_n50'], '.0f')} | {rank_display['R3']} | "
            f"{rank_display['R4']} | {rank_display['R5']} | "
            f"{row['top2_layers']} |"
        )
    counted_ranks = ", ".join(
        f"{layer} rank {fmt(errhmm[f'rank_{layer}'], '.1f')}"
        for layer in counted
    )
    lines += [
        "",
        "## Decision",
        "",
        f"GCerrHMM places top-2 in **{errhmm['top2_layers']} of {len(counted)}** "
        f"layers with a real anchor ({counted_ranks}).",
        "",
    ]
    if excluded:
        lines.insert(
            6,
            "Layers excluded from the ranking (no real anchor exists): "
            + ", ".join(sorted(excluded)),
        )
    if headline_branch:
        lines.append(
            "**Headline: downstream consistency on a GC-heterogeneous eukaryotic "
            "genome.** GCerrHMM's distances to the real anchor place it in the "
            "top-2 of at least two layers."
        )
    else:
        lines.append(
            "**Headline: compositional fidelity on GC-heterogeneous genomes plus a "
            "trainable framework.** GCerrHMM does not place top-2 in two layers. "
            "Report the R3/R5 distance ranks above as the fallback evidence."
        )
        if "R4" in excluded:
            lines.append(
                "Protocol note: the pre-registered fallback sentence naming "
                "PBSim3-sample's real error spectrum is not evaluated here because "
                "R4 is excluded; replace it with the observed R3/R5 rank table."
            )
    lines += [
        "",
        "## Provenance",
        "",
        "## Reporting commitments",
        "",
        "1. The layer where GCerrHMM is not top-2 is named in the abstract-level "
        "narrative, not hidden in the supplement.",
        "2. A split between the primary metric and the sub-metrics (for example a "
        "top-2 overall F1 driven by SNPs while the indel F1 is last) is reported "
        "as the result rather than resolved in the friendlier direction.",
    ]
    provenance = []
    for key, _ in TOOLS:
        for name in (
            f"mapping_{key}_{SPECIES}_{SUBSPECIES}.json",
            f"variant_{key}_{SPECIES}_{SUBSPECIES}.json",
            f"assembly_flye_{key}_{SPECIES}_{SUBSPECIES}.json",
        ):
            path = next(
                (directory / name for directory in dirs if (directory / name).is_file()),
                None,
            )
            if path is not None:
                provenance.append(f"- {name} sha256: `{sha256_file(path)}`")
    provenance_index = lines.index("## Provenance") + 2
    lines[provenance_index:provenance_index] = provenance + [""]

    md_path = out_dir / "delta_to_real_decision.md"
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print(f"\nwritten: {csv_path}")
    print(f"written: {md_path}")


if __name__ == "__main__":
    main()
