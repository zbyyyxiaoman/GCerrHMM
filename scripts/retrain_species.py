#!/usr/bin/env python3
"""Retrain each species at several MAPQ thresholds and tabulate the effect.

Produces the transparency table requested for the supplementary material:
for every species and threshold, how many reads survive the filter and what
error rate the resulting model implies. Also writes the trained models so the
selected threshold can be promoted and the alternatives kept for the
sensitivity analysis.

Usage:
  python3 retrain_species.py --root ~/errhmm_project \
      --species Ecoli Scerevisiae ... --thresholds 0 20 --output table.csv
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
from pathlib import Path

import numpy as np


def train(code_dir: Path, bam: Path, ref: Path, out: Path, mapq: int, nongc: bool) -> bool:
    if out.exists() and out.stat().st_size > 0:
        return True
    command = [
        "python3", str(code_dir / "src" / "train_errhmm.py"),
        "--bam", str(bam), "--ref", str(ref), "--output", str(out),
        "--gc-bins", "10", "--max-reads", "20000", "--min-mapq", str(mapq),
    ]
    if nongc:
        command.append("--nongc")
    result = subprocess.run(command, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  training failed (mapq={mapq}, nongc={nongc}): {result.stderr[-400:]}")
        return False
    return True


def summarise(model_path: Path) -> dict:
    payload = json.loads(model_path.read_text(encoding="utf-8"))
    meta = payload.get("training_metadata", {})
    counts = np.asarray(meta["global_transition_counts"], dtype=float)
    row = counts[0]
    total = row.sum() or 1.0
    return {
        "reads_total": meta.get("reads_total"),
        "reads_processed": meta.get("reads_processed"),
        "match_self": round(row[0] / total, 6),
        "substitution": round(row[1] / total, 6),
        "insertion_event": round(row[2] / total, 6),
        "deletion_event": round(row[3:].sum() / total, 6),
        "implied_event_rate": round(1 - row[0] / total, 6),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    parser.add_argument("--species", nargs="+", required=True)
    parser.add_argument("--thresholds", nargs="+", type=int, default=[0, 20])
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--model-dir",
        default=None,
        help="where to write the trained models (default: "
             "results/retrain_chain/models)",
    )
    args = parser.parse_args()

    root = Path(args.root).expanduser()
    code_dir = root / "code"
    model_dir = (
        Path(args.model_dir).expanduser()
        if args.model_dir
        else root / "results" / "retrain_chain" / "models"
    )
    model_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    for species in args.species:
        bam = root / "data" / "real_reads_verified" / f"{species}_ont_aligned.bam"
        ref = root / "data" / "references" / f"{species}_ref.fa"
        # Fail loudly. Silently skipping a species leaves its model on the old
        # threshold while everything else moves on, which produces a mixed
        # model set that no downstream table can detect.
        missing = [str(path) for path in (bam, ref) if not path.exists()]
        if missing:
            raise SystemExit(
                f"{species}: missing required input(s): {', '.join(missing)}. "
                "Align the real reads first, or remove the species from the "
                "list explicitly."
            )
        for threshold in args.thresholds:
            model = model_dir / f"{species}_errhmm_mapq{threshold}.json"
            print(f"[{species}] training min-mapq={threshold}")
            if not train(code_dir, bam, ref, model, threshold, nongc=False):
                continue
            stats = summarise(model)
            retention = (
                stats["reads_processed"] / stats["reads_total"]
                if stats["reads_total"] else None
            )
            rows.append({
                "species": species,
                "min_mapq": threshold,
                "reads_total": stats["reads_total"],
                "reads_processed": stats["reads_processed"],
                "retention": round(retention, 4) if retention else "",
                "substitution": stats["substitution"],
                "insertion_event": stats["insertion_event"],
                "deletion_event": stats["deletion_event"],
                "implied_event_rate": stats["implied_event_rate"],
            })
            # The GC-unaware control is needed downstream as well.
            nongc = model_dir / f"{species}_errhmm_1bin_mapq{threshold}.json"
            train(code_dir, bam, ref, nongc, threshold, nongc=True)

    fields = [
        "species", "min_mapq", "reads_total", "reads_processed", "retention",
        "substitution", "insertion_event", "deletion_event", "implied_event_rate",
    ]
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {out} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
