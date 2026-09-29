#!/usr/bin/env python3
"""Summarize repeated HiFi GCerrHMM cross-platform runs."""

import argparse
import csv
import json
import statistics
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--coverage", type=float, default=10.0)
    args = parser.parse_args()

    rows = []
    for path in sorted(Path(args.input_root).glob("**/gc_improvement_summary.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        metadata = payload.get("metadata", {})
        if float(metadata.get("coverage", 0)) != args.coverage:
            continue
        deltas = payload.get("gc_minus_1bin", {})
        rows.append({
            "species": metadata.get("species", ""),
            "coverage": metadata.get("coverage", ""),
            "seed": metadata.get("seed", ""),
            "composite_delta": deltas.get("composite", ""),
            "kmer_delta": deltas.get("kmer", ""),
            "gc_delta": deltas.get("gc", ""),
            "read_length_delta": deltas.get("read_length", ""),
            "qv_delta": deltas.get("qv", ""),
            "source": str(path),
        })
    if not rows:
        raise SystemExit("no HiFi seed result files found")

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    md = output.with_suffix(".md")
    lines = [
        "# HiFi cross-platform seed summary",
        "",
        "| Seed | Composite Δ | k-mer Δ | GC Δ | Read-length Δ |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['seed']} | {float(row['composite_delta']):+.2f} | "
            f"{float(row['kmer_delta']):+.2f} | {float(row['gc_delta']):+.2f} | "
            f"{float(row['read_length_delta']):+.2f} |"
        )
    for field, label in (
        ("composite_delta", "Composite Δ mean±SD"),
        ("kmer_delta", "k-mer Δ mean±SD"),
        ("gc_delta", "GC Δ mean±SD"),
    ):
        values = [float(row[field]) for row in rows]
        lines.append(
            f"\n{label}: {statistics.fmean(values):+.3f} ± "
            f"{statistics.stdev(values) if len(values) > 1 else 0.0:.3f}"
        )
    md.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
