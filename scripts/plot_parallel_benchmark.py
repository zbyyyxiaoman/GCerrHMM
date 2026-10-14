#!/usr/bin/env python3
"""Plot parallel simulator scaling from parallel_benchmark.csv."""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    with Path(args.input).open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    threads = np.asarray([int(row["threads_requested"]) for row in rows])
    seconds = np.asarray([float(row["wall_seconds"]) for row in rows])
    speedup = np.asarray([float(row["speedup_vs_1_thread"]) for row in rows])

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), constrained_layout=True)
    axes[0].plot(threads, seconds, marker="o", color="#4C78A8", linewidth=2)
    axes[0].set_xlabel("Requested worker processes")
    axes[0].set_ylabel("Wall time (s)")
    axes[0].set_title("Parallel simulation time")
    axes[0].grid(alpha=0.2)

    axes[1].plot(threads, speedup, marker="o", color="#2E8B57", linewidth=2, label="Measured")
    axes[1].plot(threads, threads, linestyle="--", color="#888888", label="Ideal")
    axes[1].set_xlabel("Requested worker processes")
    axes[1].set_ylabel("Speedup vs one worker")
    axes[1].set_title("Scaling and deterministic equivalence")
    axes[1].legend(frameon=False)
    axes[1].grid(alpha=0.2)
    fig.suptitle("GCerrHMM process-level parallel simulation", fontsize=13, fontweight="bold")
    for suffix in ("png", "pdf"):
        fig.savefig(output_dir / f"parallel_simulation_scaling.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(output_dir / "parallel_simulation_scaling.png")


if __name__ == "__main__":
    main()
