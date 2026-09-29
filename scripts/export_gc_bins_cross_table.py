#!/usr/bin/env python3
"""Export an independent multi-species GC-bin cross table for fig4."""

import argparse
import csv
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--stats-dir", required=True)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()

    stats_dir = Path(args.stats_dir)
    candidates = sorted(stats_dir.glob("profile_matched_gc_bins_seeds_*.csv"))
    if not candidates:
        raise SystemExit("No profile_matched_gc_bins_seeds_*.csv found")

    rows = list(csv.DictReader(candidates[-1].open()))
    csv_path = stats_dir / (
        f"gc_bins_cross_species_with_1bin_seeds_{args.run_id}.csv"
    )
    md_path = stats_dir / (
        f"gc_bins_cross_species_with_1bin_seeds_{args.run_id}.md"
    )
    fields = ["species", "gc_bins", "n", "mean", "std"]
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    with md_path.open("w") as handle:
        handle.write("# GC-bin cross-species table (3 seeds)\n\n")
        handle.write("| Species | GC bins | n | Mean | SD |\n")
        handle.write("|---|---:|---:|---:|---:|\n")
        for row in rows:
            handle.write(
                f"| {row['species']} | {row['gc_bins']} | {row['n']} | "
                f"{float(row['mean']):.2f} | {float(row['std']):.2f} |\n"
            )
    print(f"GC_BINS_CROSS_CSV={csv_path}")
    print(f"GC_BINS_CROSS_MD={md_path}")


if __name__ == "__main__":
    main()
