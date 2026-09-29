#!/usr/bin/env python3
"""Evaluate one simulated FASTQ with the project Level-1 composite score."""

import argparse
import json
from pathlib import Path

from evaluate_level1 import compare_real_vs_simulated


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--real", required=True)
    parser.add_argument("--sim", required=True)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--max-reads",
        type=int,
        default=None,
        help="Reservoir size for read-length/QV/GC statistics.",
    )
    parser.add_argument(
        "--kmer-reads",
        type=int,
        default=None,
        help="Maximum reads used for the k-mer spectrum.",
    )
    parser.add_argument(
        "--kmer-sample-size",
        type=int,
        default=None,
        help="Maximum k-mer vocabulary retained for correlation.",
    )
    parser.add_argument(
        "--cache-dir",
        default=None,
        help="Directory for the deterministic real-FASTQ Level-1 cache.",
    )
    args = parser.parse_args()

    result = compare_real_vs_simulated(
        args.real,
        args.sim,
        args.ref,
        max_reads=args.max_reads,
        kmer_max_reads=args.kmer_reads,
        kmer_sample_size=args.kmer_sample_size,
        cache_dir=args.cache_dir,
    )
    Path(args.output).write_text(
        json.dumps(result, indent=2),
        encoding="utf-8",
    )
    composite = result.get("composite_score", {})
    print(
        f"SIM_EVAL_DONE output={args.output} "
        f"composite_score={composite.get('composite_score', 0)} "
        f"sub_scores={composite.get('sub_scores', {})}"
    )


if __name__ == "__main__":
    main()
