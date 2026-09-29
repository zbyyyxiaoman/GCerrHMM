#!/usr/bin/env python3
"""Recompute route-level QV similarity with Wasserstein distance."""

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import wasserstein_distance


def histogram_distance(real_hist, sim_hist):
    real_counts = np.asarray(real_hist[0], dtype=float)
    sim_counts = np.asarray(sim_hist[0], dtype=float)
    edges = np.asarray(real_hist[1], dtype=float)
    centers = (edges[:-1] + edges[1:]) / 2.0
    if real_counts.sum() == 0 or sim_counts.sum() == 0:
        return 0.0, 0.0
    distance = float(wasserstein_distance(
        centers, centers,
        u_weights=real_counts / real_counts.sum(),
        v_weights=sim_counts / sim_counts.sum(),
    ))
    score = 100.0 * math.exp(-distance / 10.0)
    return distance, score


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tables-dir", required=True)
    parser.add_argument("--stats-dir", required=True)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()

    tables = Path(args.tables_dir)
    stats = Path(args.stats_dir)
    rows = []
    grouped = defaultdict(list)
    for path in sorted(tables.glob("level1_*.json")):
        key = path.stem.removeprefix("level1_")
        if "_route_" not in key or "_r" not in key:
            continue
        species, route_and_rep = key.split("_route_", 1)
        route_base, rep_text = route_and_rep.rsplit("_r", 1)
        if not rep_text.isdigit():
            continue
        data = json.loads(path.read_text())
        distance, score = histogram_distance(
            data["real"]["qv"]["histogram"],
            data["simulated"]["qv"]["histogram"],
        )
        row = {
            "species": species,
            "route": "route_" + route_base,
            "replicate": int(rep_text),
            "wasserstein_distance": distance,
            "wasserstein_score": score,
        }
        rows.append(row)
        grouped[(species, row["route"])].append(score)

    detail_csv = stats / f"qv_wasserstein_{args.run_id}.csv"
    summary_csv = stats / f"qv_wasserstein_summary_{args.run_id}.csv"
    summary_md = stats / f"qv_wasserstein_summary_{args.run_id}.md"
    with detail_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["species", "route", "replicate", "wasserstein_distance", "wasserstein_score"],
        )
        writer.writeheader()
        writer.writerows(rows)

    summary_rows = []
    for (species, route), values in sorted(grouped.items()):
        summary_rows.append({
            "species": species,
            "route": route,
            "n": len(values),
            "mean_score": float(np.mean(values)),
            "sd_score": float(np.std(values, ddof=1)) if len(values) > 1 else 0.0,
        })
    with summary_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["species", "route", "n", "mean_score", "sd_score"]
        )
        writer.writeheader()
        writer.writerows(summary_rows)
    with summary_md.open("w") as handle:
        handle.write("# QV Wasserstein summary\n\n")
        handle.write("| Species | Route | n | Mean score | SD |\n")
        handle.write("|---|---|---:|---:|---:|\n")
        for row in summary_rows:
            handle.write(
                f"| {row['species']} | {row['route']} | {row['n']} | "
                f"{row['mean_score']:.2f} | {row['sd_score']:.2f} |\n"
            )
    print(f"QV_WASSERSTEIN_DETAIL={detail_csv}")
    print(f"QV_WASSERSTEIN_SUMMARY={summary_md}")


if __name__ == "__main__":
    main()
