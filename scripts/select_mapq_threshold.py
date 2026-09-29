#!/usr/bin/env python3
"""Pick the training MAPQ threshold by error-spectrum distance.

Criterion, fixed before the numbers were seen: choose the model whose
simulated (substitution, insertion, deletion) rate vector is closest to the
real read set in Euclidean distance, after normalising each axis by the real
rate so that no single error class dominates. Ties are broken by the
GC-stratified curve distance.

Usage:
  python3 select_mapq_threshold.py --input-dir results/mapq_threshold_selection \
      --output selection.json
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

AXES = ("substitution_rate", "insertion_rate", "deletion_rate")


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def normalised_distance(real: dict, sim: dict) -> float:
    total = 0.0
    for axis in AXES:
        base = real[axis] or 1e-9
        total += ((sim[axis] - real[axis]) / base) ** 2
    return math.sqrt(total)


def gc_curve_distance(real: dict, sim: dict) -> float:
    """Mean absolute difference over GC bins that are populated on both sides."""
    total = 0.0
    bins = 0
    for real_bin, sim_bin in zip(real["by_gc"], sim["by_gc"]):
        if real_bin["error_rate"] is None or sim_bin["error_rate"] is None:
            continue
        if real_bin["aligned_bases"] < 1000 or sim_bin["aligned_bases"] < 1000:
            continue
        total += abs(sim_bin["error_rate"] - real_bin["error_rate"])
        bins += 1
    return total / bins if bins else float("nan")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    root = Path(args.input_dir)
    real = load(root / "spectrum_real.json")

    candidates = []
    for path in sorted(root.glob("spectrum_mapq*.json")):
        sim = load(path)
        candidates.append(
            {
                "model": path.stem.replace("spectrum_", ""),
                "file": path.name,
                "global": sim["global"],
                "reads": sim["reads"],
                "normalised_distance": round(
                    normalised_distance(real["global"], sim["global"]), 4
                ),
                "gc_curve_distance": round(gc_curve_distance(real, sim), 4),
            }
        )

    candidates.sort(key=lambda row: (row["normalised_distance"], row["gc_curve_distance"]))
    payload = {
        "criterion": (
            "minimise normalised Euclidean distance in "
            "(substitution, insertion, deletion) rate space; "
            "ties broken by mean absolute GC-stratified error difference"
        ),
        "reference": {"label": real["label"], "global": real["global"]},
        "candidates": candidates,
        "selected": candidates[0]["model"] if candidates else None,
    }

    out = Path(args.output)
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
