#!/usr/bin/env python3
"""GC-conditional error fidelity with the simulated side pooled.

The 50-point curve correlation is attenuated by sampling noise on the simulated
side: each GC cell of a single 10x replicate holds few aligned bases, so the
correlation between a precise real curve and a noisy simulated one is biased
towards zero (errors-in-variables attenuation).

Pooling the replicates raises the simulated coverage without touching the
estimator, the GC cells or the decision rule - it is variance reduction, the
same category as raising the k-mer sample. Replicates are summed before any
rate is formed, so the result is one curve estimated from 20x (panel) or 30x
(three ablation seeds).

usage:
  python3 pooled_gc_curve.py --ref ref.fa --real-bam real.bam \
      --sim-bam run1.bam --sim-bam run2.bam --curve-bins 50 --output out.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from gc_error_fidelity import (  # noqa: E402
    accumulate_by_window,
    longest_contig,
    summarise_curve,
    summarise_windows,
    window_gc,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", required=True)
    parser.add_argument("--real-bam", required=True)
    parser.add_argument("--sim-bam", action="append", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--contig", default=None)
    parser.add_argument("--window-size", type=int, default=100)
    parser.add_argument("--curve-bins", type=int, default=50)
    parser.add_argument("--curve-min-bases", type=int, default=2000)
    parser.add_argument("--bootstrap", type=int, default=200)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    ref = Path(args.ref)
    contig = args.contig or longest_contig(ref)
    real = accumulate_by_window(Path(args.real_bam), ref, contig,
                                args.window_size)
    pooled = None
    for bam in args.sim_bam:
        counts = accumulate_by_window(Path(bam), ref, contig, args.window_size)
        if pooled is None:
            pooled = counts
        else:
            pooled = {
                "aligned": pooled["aligned"] + counts["aligned"],
                "errors": pooled["errors"] + counts["errors"],
            }
    gc_frac = window_gc(ref, contig, args.window_size)
    curve = summarise_curve(real, pooled, gc_frac, curve_bins=args.curve_bins,
                            min_bases=args.curve_min_bases,
                            bootstrap=args.bootstrap)
    windows = summarise_windows(real, pooled, min_support=10)
    payload = {
        "label": args.label,
        "replicates_pooled": len(args.sim_bam),
        "sim_bams": args.sim_bam,
        "contig": contig,
        "curve": curve,
        "window_level": windows,
    }
    Path(args.output).write_text(json.dumps(payload, indent=2))
    print(json.dumps({k: v for k, v in curve.items() if k != "per_point"},
                     indent=2))


if __name__ == "__main__":
    main()
