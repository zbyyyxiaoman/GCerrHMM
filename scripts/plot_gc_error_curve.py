#!/usr/bin/env python3
"""Draw the GC-error curve (Fig A) and annotate every dataset with its fidelity.

The curve is the primary evidence for the GC-conditioning claim, and a curve
without a number is an illustration. For each dataset the script reports the
summaries that `docs/gc_claim_evidence_chain.md` pre-registers - Pearson r,
Spearman rho and mean absolute deviation against the real per-GC-bin error
rate - through the shared implementation in `scripts/gc_error_fidelity.py`, so
the figure and the table cannot drift apart.

usage:
  python3 plot_gc_error_curve.py --ref ref.fa --real-bam real.bam \
      --dataset errHMM=sim.bam --dataset NanoSim=nano.bam \
      --output-prefix docs/framework_figures/figA_gc_error_ecoli
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pysam  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gc_error_fidelity import (  # noqa: E402
    accumulate, accumulate_by_window, summarise, summarise_curve,
    summarise_windows, window_gc,
)

COLOURS = ["#d1495b", "#00798c", "#edae49", "#30638e", "#6a994e", "#8d6a9f"]


def longest_contig(ref: Path) -> str:
    with pysam.FastaFile(str(ref)) as fa:
        return max(zip(fa.references, fa.lengths), key=lambda x: x[1])[0]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", required=True)
    parser.add_argument("--real-bam", required=True)
    parser.add_argument("--dataset", action="append", required=True,
                        help="LABEL=BAM (repeatable)")
    parser.add_argument("--contig", default=None)
    parser.add_argument("--bins", type=int, default=10)
    parser.add_argument("--window-size", type=int, default=100)
    parser.add_argument("--min-support", type=int, default=20000)
    parser.add_argument("--curve-bins", type=int, default=50,
                        help="GC cells for the higher-power curve correlation")
    parser.add_argument("--output-prefix", required=True)
    args = parser.parse_args()

    ref = Path(args.ref)
    contig = args.contig or longest_contig(ref)
    real = accumulate(Path(args.real_bam), ref, contig, args.window_size,
                      args.bins)
    real_summary = summarise(real, real, args.min_support)
    real_win = accumulate_by_window(Path(args.real_bam), ref, contig,
                                    args.window_size)
    gc_frac = window_gc(ref, contig, args.window_size)

    curves = []
    for spec in args.dataset:
        label, _, bam = spec.partition("=")
        if not bam:
            raise SystemExit(f"bad --dataset {spec!r}, expected LABEL=BAM")
        sim = accumulate(Path(bam), ref, contig, args.window_size, args.bins)
        bin_summary = summarise(real, sim, args.min_support)
        sim_win = accumulate_by_window(Path(bam), ref, contig, args.window_size)
        win_summary = summarise_windows(real_win, sim_win, min_support=10)
        curve_summary = summarise_curve(real_win, sim_win, gc_frac,
                                        curve_bins=args.curve_bins)
        curves.append((label, bin_summary, win_summary, curve_summary))

    prefix = Path(args.output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    with prefix.with_suffix(".csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["dataset", "gc_bin", "aligned_bases", "error_rate",
                         "pearson_r", "spearman_rho", "mad", "bins_used",
                         "window_r", "window_ci_low", "window_ci_high",
                         "windows_used"])
        for label, summary in [("real", real_summary),
                               *[(c[0], c[1]) for c in curves]]:
            for row in summary["per_bin"]:
                writer.writerow([label, row["gc_bin"], row["sim_aligned"],
                                 row["sim_error_rate"], summary["pearson_r"],
                                 summary["spearman_rho"], summary["mad"],
                                 summary.get("bins_used"), "", "", "", ""])
        for label, _bin, win, cur in curves:
            writer.writerow([label, "window_level", "", "", "", "", win["mad"],
                             "", win["pearson_r"], win["ci_low"],
                             win["ci_high"], win["windows_used"]])
            writer.writerow([label, "curve_level", "", "", "", "", cur["mad"],
                             "", cur["pearson_r"], cur["ci_low"],
                             cur["ci_high"], cur["points_used"]])

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    used = [row for row in real_summary["per_bin"] if row["used"]]
    ax.plot([r["gc_bin"] for r in used],
            [r["real_error_rate"] * 100 for r in used],
            color="black", linewidth=2.6, marker="o", label="real ONT",
            zorder=5)

    for idx, (label, summary, win, cur) in enumerate(curves):
        rows = [r for r in summary["per_bin"]
                if r["used"] and r["sim_error_rate"] is not None]
        if not rows:
            continue
        r = summary["pearson_r"]
        if cur["pearson_r"] is not None:
            tag = (f"50-pt r={cur['pearson_r']:.2f} "
                   f"[{cur['ci_low']:.2f},{cur['ci_high']:.2f}]")
        elif r is not None:
            tag = f"bin r={r:.2f} (n={summary.get('bins_used', 0)})"
        else:
            tag = f"r=n/a, MAD={summary['mad'] * 100:.2f}pp"
        ax.plot([row["gc_bin"] for row in rows],
                [row["sim_error_rate"] * 100 for row in rows],
                color=COLOURS[idx % len(COLOURS)], linewidth=1.8,
                marker="s", markersize=4, label=f"{label} ({tag})")

    ax.set_xlabel("reference GC bin (10 % widths)")
    ax.set_ylabel("error rate (%)")
    ax.set_title("Error rate by local GC content")
    ax.grid(alpha=0.25, linewidth=0.5)
    ax.set_axisbelow(True)
    ax.legend(fontsize=7.5, frameon=False)
    fig.tight_layout()
    for suffix in ("png", "pdf"):
        fig.savefig(prefix.with_suffix(f".{suffix}"), dpi=300)
    plt.close(fig)

    print(f"{'dataset':<22} {'8-bin r':>7} {'n':>3} {'50-pt r':>8} "
          f"{'50-pt CI':>18} {'n_pt':>5} {'win r':>7} {'MAD(pp)':>8}")
    for label, summary, win, cur in curves:
        r = summary["pearson_r"]
        cr = cur["pearson_r"]
        cci = ("n/a" if cr is None
               else f"[{cur['ci_low']:.3f},{cur['ci_high']:.3f}]")
        wr = win.get("pearson_r")
        print(f"{label:<22} "
              f"{'n/a' if r is None else format(r, '.3f'):>7} "
              f"{summary.get('bins_used', 0):>3} "
              f"{'n/a' if cr is None else format(cr, '.3f'):>8} "
              f"{cci:>18} {cur.get('points_used', 0):>5} "
              f"{'n/a' if wr is None else format(wr, '.3f'):>7} "
              f"{summary['mad'] * 100:>8.3f}")
    print(f"written: {prefix.with_suffix('.png')}")


if __name__ == "__main__":
    main()
