#!/usr/bin/env python3
"""Benchmark the process-parallel GCerrHMM read generator."""

import argparse
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path


def run(command, log):
    log.write("$ " + " ".join(str(item) for item in command) + "\n")
    log.flush()
    result = subprocess.run(
        [str(item) for item in command],
        stdout=log,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"command failed: {command}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-dir", default=".")
    parser.add_argument("--repo-dir", default=None)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--coverage", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--threads", "--workers",
        dest="threads",
        default="1,2,4,8",
        help="Comma-separated process-worker counts.",
    )
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    root = Path(args.project_dir).resolve()
    repo = Path(args.repo_dir).resolve() if args.repo_dir else root
    output_dir = Path(args.output_dir).resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise SystemExit(f"refusing to overwrite: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "parallel_benchmark.log"
    thread_counts = [int(item) for item in args.threads.split(",") if item.strip()]
    rows = []

    with log_path.open("w", encoding="utf-8") as log:
        for thread_count in thread_counts:
            prefix = output_dir / f"threads_{thread_count}"
            command = [
                sys.executable, repo / "src/generate_simulated.py",
                "--ref", Path(args.ref).resolve(),
                "--output-prefix", prefix,
                "--coverage", str(args.coverage),
                "--platform", "ont",
                "--errhmm-model", Path(args.model).resolve(),
                "--seed", str(args.seed),
                "--workers", str(thread_count),
            ]
            if shutil.which("nice"):
                command = ["nice", "-n", "19"] + command
            run(command, log)
            stats = json.loads(Path(f"{prefix}.stats.json").read_text(encoding="utf-8"))
            rows.append({
                "threads_requested": thread_count,
                "threads_used": stats.get("threads_used"),
                "wall_seconds": stats.get("wall_seconds"),
                "reads_per_second": stats.get("reads_per_second"),
                "num_reads": stats.get("num_reads"),
                "total_bases": stats.get("total_bases"),
                "mean_length": stats.get("mean_length"),
                "n50": stats.get("n50"),
            })

    baseline = rows[0]["wall_seconds"] if rows else 0.0
    for row in rows:
        row["speedup_vs_1_thread"] = (
            baseline / row["wall_seconds"] if row["wall_seconds"] else 0.0
        )

    csv_path = output_dir / "parallel_benchmark.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    md_path = output_dir / "parallel_benchmark.md"
    lines = [
        "# Parallel simulation benchmark",
        "",
        "| Threads requested | Threads used | Wall seconds | Reads/s | Speedup vs 1 |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['threads_requested']} | {row['threads_used']} | "
            f"{row['wall_seconds']:.3f} | {row['reads_per_second']:.1f} | "
            f"{row['speedup_vs_1_thread']:.2f}x |"
        )
    lines += [
        "",
        "All runs use the same reference, model, coverage, and seed. Read-count,",
        "total-base, mean-length, and N50 fields are retained for equivalence",
        "checking; the generator derives per-read seeds so worker count does not",
        "change those distribution statistics.",
    ]
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(csv_path)


if __name__ == "__main__":
    main()
