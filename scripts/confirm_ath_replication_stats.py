#!/usr/bin/env python3
"""Summarise the pre-specified second A. thaliana replication cohort."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy import stats


PHASE1_SEEDS = (44, 55, 66)
PHASE2_SEEDS = (77, 88, 99, 111, 122, 133)
OLD_SEEDS = (11, 22, 33)


def load_r(path: Path) -> float:
    payload = json.loads(path.read_text())
    value = payload.get("pearson_r")
    if value is None:
        value = payload.get("curve", {}).get("pearson_r")
    if value is None:
        raise ValueError(f"no Pearson r in {path}")
    return float(value)


def load_status(status_dir: Path, seeds: tuple[int, ...]) -> dict[int, dict[int, float]]:
    values = {1: {}, 5: {}}
    for seed in seeds:
        for bins in (1, 5):
            path = status_dir / f"ablation_Athaliana_{bins}bin_seed{seed}.gc.json"
            values[bins][seed] = load_r(path)
    return values


def load_old(summary: Path) -> dict[int, dict[int, float]]:
    values = {1: {}, 5: {}}
    with summary.open(newline="") as handle:
        for row in csv.DictReader(handle):
            for bins in (1, 5):
                for seed in OLD_SEEDS:
                    if row["label"] == f"ablation_Athaliana_{bins}bin_seed{seed}":
                        values[bins][seed] = float(row["gc_pearson_r"])
    return values


def require_complete(
    values: dict[int, dict[int, float]],
    seeds: tuple[int, ...],
    label: str,
) -> None:
    missing = [
        f"{bins}bin:{seed}"
        for bins in (1, 5)
        for seed in seeds
        if seed not in values[bins]
    ]
    if missing:
        raise ValueError(f"{label} is missing required values: {', '.join(missing)}")


def paired_summary(values: dict[int, dict[int, float]], seeds: tuple[int, ...]) -> dict:
    y = np.asarray([values[5][seed] for seed in seeds], dtype=float)
    x = np.asarray([values[1][seed] for seed in seeds], dtype=float)
    delta = y - x
    n = len(seeds)
    mean = float(delta.mean())
    sd = float(delta.std(ddof=1)) if n > 1 else 0.0
    se = sd / np.sqrt(n) if n > 1 else 0.0
    if n > 1:
        critical = float(stats.t.ppf(0.975, df=n - 1))
        test = stats.ttest_rel(y, x)
        low, high = mean - critical * se, mean + critical * se
        t_stat, p_value = float(test.statistic), float(test.pvalue)
    else:
        low = high = None
        t_stat = p_value = None
    wilcoxon = (
        float(stats.wilcoxon(delta, alternative="two-sided").pvalue)
        if np.any(delta != 0) else None
    )
    return {
        "seeds": list(seeds),
        "n": n,
        "r_1bin": {str(seed): float(values[1][seed]) for seed in seeds},
        "r_5bin": {str(seed): float(values[5][seed]) for seed in seeds},
        "delta": {str(seed): float(value) for seed, value in zip(seeds, delta)},
        "mean_delta": mean,
        "sd_delta": sd,
        "ci95_low": None if low is None else float(low),
        "ci95_high": None if high is None else float(high),
        "paired_t": t_stat,
        "paired_t_p": p_value,
        "wilcoxon_p": wilcoxon,
        "supported": bool(mean > 0 and p_value is not None and p_value < 0.05),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase1-status", type=Path, required=True)
    parser.add_argument("--phase2-status", type=Path, required=True)
    parser.add_argument("--old-summary", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    args = parser.parse_args()

    phase1 = load_status(args.phase1_status, PHASE1_SEEDS)
    phase2 = load_status(args.phase2_status, PHASE2_SEEDS)
    old = load_old(args.old_summary)
    require_complete(phase1, PHASE1_SEEDS, "phase-1 cohort")
    require_complete(phase2, PHASE2_SEEDS, "phase-2 cohort")
    require_complete(old, OLD_SEEDS, "discovery cohort")
    all_values = {
        1: {**old[1], **phase1[1], **phase2[1]},
        5: {**old[5], **phase1[5], **phase2[5]},
    }

    result = {
        "second_cohort": paired_summary(phase2, PHASE2_SEEDS),
        "first_cohort": paired_summary(phase1, PHASE1_SEEDS),
        "discovery_cohort": paired_summary(old, OLD_SEEDS),
        "all_seed_sensitivity": paired_summary(
            all_values,
            OLD_SEEDS + PHASE1_SEEDS + PHASE2_SEEDS,
        ),
        "provenance": {
            "phase1_status": str(args.phase1_status),
            "phase2_status": str(args.phase2_status),
            "old_summary": str(args.old_summary),
            "primary_endpoint": "paired 5-bin minus 1-bin GC-curve Pearson r",
        },
    }

    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, indent=2) + "\n")

    second = result["second_cohort"]
    first = result["first_cohort"]
    pooled = result["all_seed_sensitivity"]
    lines = [
        "# Second pre-specified A. thaliana replication cohort",
        "",
        "Primary endpoint: paired `5-bin - 1-bin` Pearson r.",
        "",
        "Protocol: `docs/confirmatory_ath_bin_replication_phase2_protocol.md`.",
        "",
        "## Second cohort (seeds 77, 88, 99, 111, 122, 133)",
        "",
        f"- seed-wise deltas: `{second['delta']}`",
        f"- mean delta: `{second['mean_delta']:.10f}`",
        f"- sample SD: `{second['sd_delta']:.10f}`",
        f"- 95% t CI: `[{second['ci95_low']:.10f}, {second['ci95_high']:.10f}]`",
        f"- paired t: `{second['paired_t']:.10f}`, "
        f"p = `{second['paired_t_p']:.10f}`",
        f"- Wilcoxon p: `{second['wilcoxon_p']:.10f}`",
        f"- pre-specified support: `{second['supported']}`",
        "",
        "| Seed | 1-bin r | 5-bin r | Delta |",
        "|---:|---:|---:|---:|",
    ]
    for seed in second["seeds"]:
        key = str(seed)
        lines.append(
            f"| {seed} | {second['r_1bin'][key]:.10f} | "
            f"{second['r_5bin'][key]:.10f} | "
            f"{second['delta'][key]:.10f} |"
        )
    lines += [
        "",
        "## First cohort (seeds 44, 55, 66)",
        "",
        f"- mean delta: `{first['mean_delta']:.10f}`",
        f"- paired t p: `{first['paired_t_p']:.10f}`",
        f"- support: `{first['supported']}`",
        "",
        "## Pooled 12-seed sensitivity",
        "",
        f"- mean delta: `{pooled['mean_delta']:.10f}`",
        f"- 95% t CI: `[{pooled['ci95_low']:.10f}, "
        f"{pooled['ci95_high']:.10f}]`",
        f"- paired t p: `{pooled['paired_t_p']:.10f}`",
        "",
        "## Endpoint implementation",
        "",
        "The endpoint-only implementation uses the same reference windows, "
        "10 equal-width GC bins, CIGAR event definition and minimum-support "
        "rule as the complete fidelity script. Before the second cohort was "
        "scored, it reproduced the six first-cohort top-level `pearson_r` "
        "values exactly; all six absolute differences were `0.0`. The unused "
        "window-level bootstrap and descriptive curve were omitted from the "
        "confirmation run and remain available in the complete fidelity "
        "script. The equivalence table is retained in "
        "`docs/confirmatory_ath_bin_replication_phase2_bundle/"
        "EQUIVALENCE.md`.",
        "",
        "The endpoint JSON files and their SHA-256 manifest are stored under "
        "`docs/confirmatory_ath_bin_replication_phase2_bundle/results/"
        "confirm_ath_bins_phase2_20261014/status/`.",
        "",
        "The original cross-species gate remains unchanged.",
        "",
    ]
    args.output_md.write_text("\n".join(lines))
    print(args.output_md.read_text())


if __name__ == "__main__":
    main()
