#!/usr/bin/env python3
"""Measure training and simulation wall time, peak RSS and process scaling."""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


def run_timed(command: list[str], log_path: Path, time_path: Path) -> None:
    with log_path.open("w", encoding="utf-8") as log:
        log.write("$ " + " ".join(command) + "\n")
        log.flush()
        result = subprocess.run(
            [
                "/usr/bin/time",
                "-v",
                "-o",
                str(time_path),
                *command,
            ],
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )
    if result.returncode != 0:
        raise RuntimeError(f"command failed: {' '.join(command)}")


def parse_time(path: Path) -> tuple[float, float]:
    text = path.read_text(encoding="utf-8", errors="replace")
    wall = re.search(r"Elapsed \(wall clock\) time.*?\): ([0-9:.]+)", text)
    rss = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    if not wall or not rss:
        raise RuntimeError(f"cannot parse GNU time output: {path}")
    parts = [float(part) for part in wall.group(1).split(":")]
    if len(parts) == 2:
        hours = 0.0
        minutes, seconds = parts
    elif len(parts) == 3:
        hours, minutes, seconds = parts
    else:
        raise RuntimeError(f"unexpected GNU time format in {path}: {wall.group(1)}")
    return hours * 3600 + minutes * 60 + seconds, int(rss.group(1)) / 1024


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-dir", default=".")
    parser.add_argument("--repo-dir", default=None)
    parser.add_argument("--species", default="Ecoli")
    parser.add_argument("--coverage", type=float, default=10.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--workers", default="1,2,4")
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--max-reads", type=int, default=20000)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    root = Path(args.project_dir).resolve()
    repo = Path(args.repo_dir).resolve() if args.repo_dir else root
    out = Path(args.output_dir).resolve()
    if out.exists() and any(out.iterdir()):
        raise SystemExit(f"refusing to overwrite non-empty output: {out}")
    out.mkdir(parents=True, exist_ok=True)
    logs = out / "logs"
    logs.mkdir()

    if not Path("/usr/bin/time").exists():
        raise SystemExit("GNU /usr/bin/time is required for peak-RSS measurement")

    ref = root / "data/references" / f"{args.species}_ref.fa"
    bam = root / "data/real_reads_verified" / f"{args.species}_ont_aligned.bam"
    profile = root / "data/tmp/read_profiles" / f"{args.species}.json"
    for path in (ref, bam, profile):
        if not path.exists():
            raise SystemExit(f"missing input: {path}")

    nice = ["nice", "-n", "19"] if shutil.which("nice") else []
    model = out / f"{args.species}_model_gc10.json"
    train_command = nice + [
        sys.executable,
        str(repo / "src/train_errhmm.py"),
        "--bam",
        str(bam),
        "--ref",
        str(ref),
        "--output",
        str(model),
        "--gc-bins",
        "10",
        "--window-size",
        "100",
        "--min-mapq",
        str(args.min_mapq),
        "--max-reads",
        str(args.max_reads),
        "--min-reads",
        "1000",
    ]
    run_timed(train_command, logs / "training.log", logs / "training.time")

    worker_counts = [int(item) for item in args.workers.split(",") if item.strip()]
    rows: list[dict[str, object]] = []
    for workers in worker_counts:
        prefix = out / f"sim_{workers}"
        command = nice + [
            sys.executable,
            str(repo / "src/generate_simulated.py"),
            "--ref",
            str(ref),
            "--output-prefix",
            str(prefix),
            "--coverage",
            str(args.coverage),
            "--platform",
            "ont",
            "--errhmm-model",
            str(model),
            "--profile-json",
            str(profile),
            "--require-profile",
            "--seed",
            str(args.seed),
            "--workers",
            str(workers),
        ]
        run_timed(command, logs / f"sim_{workers}.log", logs / f"sim_{workers}.time")
        wall, rss = parse_time(logs / f"sim_{workers}.time")
        stats = json.loads(Path(f"{prefix}.stats.json").read_text(encoding="utf-8"))
        rows.append({
            "stage": "simulate",
            "workers": workers,
            "wall_seconds": wall,
            "peak_rss_mb": rss,
            "reads": stats["num_reads"],
            "reads_per_second": stats["num_reads"] / wall,
            "speedup_vs_1": 1.0,
        })

    train_wall, train_rss = parse_time(logs / "training.time")
    rows.insert(0, {
        "stage": "train_gc10",
        "workers": 1,
        "wall_seconds": train_wall,
        "peak_rss_mb": train_rss,
        "reads": "",
        "reads_per_second": "",
        "speedup_vs_1": 1.0,
    })
    baseline = next(row["wall_seconds"] for row in rows if row["stage"] == "simulate")
    for row in rows:
        if row["stage"] == "simulate":
            row["speedup_vs_1"] = baseline / row["wall_seconds"]

    csv_path = out / "resource_benchmark.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "# GCerrHMM resource benchmark",
        "",
        f"Species: {args.species}; coverage: {args.coverage}; seed: {args.seed}",
        f"Training MAPQ threshold: {args.min_mapq}; workers: {args.workers}",
        "",
        "| Stage | Workers | Wall seconds | Peak RSS MB | Reads | Reads/s | Speedup |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        reads = row["reads"] if row["reads"] != "" else ""
        rps = f"{row['reads_per_second']:.1f}" if row["reads_per_second"] != "" else ""
        lines.append(
            f"| {row['stage']} | {row['workers']} | {row['wall_seconds']:.2f} | "
            f"{row['peak_rss_mb']:.1f} | {reads} | {rps} | "
            f"{row['speedup_vs_1']:.2f}x |"
        )
    md_path = out / "resource_benchmark.md"
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(md_path)


if __name__ == "__main__":
    main()
