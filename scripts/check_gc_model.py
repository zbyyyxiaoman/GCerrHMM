#!/usr/bin/env python3
"""Check that a trained model actually contains GC-conditioned structure."""

import argparse
import json
from pathlib import Path

import numpy as np


def matrix_distance(left, right):
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    return float(np.mean(np.abs(left - right)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    model = json.loads(Path(args.model).read_text(encoding="utf-8"))
    matrices = model.get("transition_probs", {})
    bins = sorted(matrices, key=lambda item: int(item))
    distances = []
    for left, right in zip(bins, bins[1:]):
        distances.append(matrix_distance(matrices[left], matrices[right]))
    mean_distance = float(np.mean(distances)) if distances else 0.0
    max_distance = float(np.max(distances)) if distances else 0.0
    nonempty = [
        key for key, matrix in matrices.items()
        if float(np.asarray(matrix, dtype=float).sum()) > 0
    ]
    payload = {
        "model": str(Path(args.model).resolve()),
        "gc_bins": int(model.get("gc_bins", len(bins))),
        "bins_with_values": len(nonempty),
        "mean_adjacent_transition_l1": mean_distance,
        "max_adjacent_transition_l1": max_distance,
        "gc_structure_present": len(bins) > 1 and max_distance > 0.0,
    }
    Path(args.output).write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    if args.strict and not payload["gc_structure_present"]:
        raise SystemExit("GC-conditioned structure was not detected")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
