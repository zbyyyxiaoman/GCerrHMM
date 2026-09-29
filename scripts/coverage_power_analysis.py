#!/usr/bin/env python3
"""How much coverage does the GC-fidelity instrument need to see an effect?

The replicated curves show a consistent direction (GC-aware above 1-bin in
every species tested) that never clears the bootstrap interval. Turning that
null result into a usable statement requires the noise floor as a function of
coverage: at 10x, 20x and 30x, how wide is the interval on the curve
correlation?

Each level is given as one or more simulated BAMs which are summed before any
rate is formed (so two 10x runs enter as 20x), and every level goes through the
same 50-point curve estimator and bootstrap as the paper's primary instrument.

usage:
  python3 coverage_power_analysis.py --ref ref.fa --real-bam real.bam \
      --level 10x=a.bam --level 20x=a.bam,b.bam --level 30x=a.bam,b.bam,c.bam \
      --label Ecoli_GCerrHMM --output results/stats/coverage_power_Ecoli.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gc_error_fidelity import (  # noqa: E402
    accumulate_by_window, longest_contig, summarise_curve, window_gc,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", required=True)
    parser.add_argument("--real-bam", required=True)
    parser.add_argument("--level", action="append", required=True,
                        help="LABEL=BAM[,BAM...] (replicates are summed)")
    parser.add_argument("--label", required=True)
    parser.add_argument("--curve-bins", type=int, default=50)
    parser.add_argument("--bootstrap", type=int, default=200)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    ref = Path(args.ref)
    contig = longest_contig(ref)
    real = accumulate_by_window(Path(args.real_bam), ref, contig, 100)
    gc_frac = window_gc(ref, contig, 100)

    rows = []
    for spec in args.level:
        label, _, bams = spec.partition("=")
        pooled = None
        for bam in bams.split(","):
            counts = accumulate_by_window(Path(bam), ref, contig, 100)
            pooled = counts if pooled is None else {
                "aligned": pooled["aligned"] + counts["aligned"],
                "errors": pooled["errors"] + counts["errors"],
            }
        curve = summarise_curve(real, pooled, gc_frac,
                                curve_bins=args.curve_bins,
                                bootstrap=args.bootstrap)
        rows.append({
            "label": args.label,
            "coverage": label,
            "replicates": len(bams.split(",")),
            "pearson_r": curve["pearson_r"],
            "ci_low": curve["ci_low"],
            "ci_high": curve["ci_high"],
            "ci_width": (None if curve["ci_low"] is None
                         else curve["ci_high"] - curve["ci_low"]),
            "points_used": curve["points_used"],
            "mad": curve["mad"],
            "real_range": curve.get("real_range"),
            "sim_range": curve.get("sim_range"),
            "real_rate_vs_gc": curve.get("real_rate_correlation_with_gc"),
        })

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0].keys())
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    xs = list(range(len(rows)))
    widths = [r["ci_width"] if r["ci_width"] is not None else float("nan")
              for r in rows]
    rs = [r["pearson_r"] if r["pearson_r"] is not None else float("nan")
          for r in rows]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(xs, widths, marker="o", color="#30638e", label="95% CI width")
    ax.plot(xs, rs, marker="s", color="#d1495b", label="curve correlation r")
    ax.set_xticks(xs)
    ax.set_xticklabels([r["coverage"] for r in rows])
    ax.set_xlabel("simulated coverage (replicates pooled)")
    ax.set_ylabel("value")
    ax.set_title(f"GC-fidelity noise floor vs coverage - {args.label}")
    ax.grid(alpha=0.25, linewidth=0.5)
    ax.set_axisbelow(True)
    ax.legend(fontsize=8, frameon=False)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(out.with_suffix(f".{suffix}"), dpi=300)
    plt.close(fig)

    print(json.dumps(rows, indent=2))
    print(f"written: {out}")


if __name__ == "__main__":
    main()
