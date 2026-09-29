#!/usr/bin/env python3
"""Build reproducible read-length and QV profiles from real FASTQ data."""

import argparse
import gzip
import json
import random
import sys
from pathlib import Path

import numpy as np
from Bio import SeqIO


def open_fastq(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--accession", required=True)
    parser.add_argument("--max-reads", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260913)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    reservoir = []
    qv_hist = np.zeros(94, dtype=np.int64)
    total_reads = 0
    for index, record in enumerate(SeqIO.parse(open_fastq(args.fastq), "fastq")):
        total_reads += 1
        if index < args.max_reads:
            reservoir.append(record)
        else:
            slot = rng.randrange(index + 1)
            if slot < args.max_reads:
                reservoir[slot] = record

    if not reservoir:
        raise SystemExit(f"No reads found in {args.fastq}")

    lengths = np.asarray([len(record.seq) for record in reservoir], dtype=float)
    for record in reservoir:
        qvs = record.letter_annotations.get("phred_quality")
        if qvs:
            arr = np.asarray(qvs, dtype=np.int64)
            qv_hist += np.bincount(arr, minlength=94)[:94]

    qv_total = int(qv_hist.sum())
    if qv_total:
        support = np.arange(94, dtype=float)
        qv_mean = float(np.average(support, weights=qv_hist))
        qv_std = float(np.sqrt(np.average((support - qv_mean) ** 2, weights=qv_hist)))
    else:
        qv_mean = qv_std = 0.0

    quantiles = np.linspace(0.0, 1.0, 101)
    profile = {
        "species": args.species,
        "accession": args.accession,
        "source_fastq": str(args.fastq),
        "source_reads": total_reads,
        "sampled_reads": len(reservoir),
        "sampling": "deterministic_reservoir",
        "read_length": {
            "mean": float(np.mean(lengths)),
            "median": float(np.median(lengths)),
            "std": float(np.std(lengths)),
            "min": int(np.min(lengths)),
            "max": int(np.max(lengths)),
            "n50": int(np.sort(lengths)[::-1][
                np.searchsorted(np.cumsum(np.sort(lengths)[::-1]), np.sum(lengths) / 2)
            ]),
            "quantile_probabilities": quantiles.tolist(),
            "quantiles": np.quantile(lengths, quantiles).tolist(),
        },
        "qv": {
            "mean": qv_mean,
            "std": qv_std,
            "histogram": qv_hist.tolist(),
        },
    }

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", encoding="utf-8") as handle:
        json.dump(profile, handle, indent=2)
    print(json.dumps({
        "species": args.species,
        "reads": total_reads,
        "sampled": len(reservoir),
        "mean_length": profile["read_length"]["mean"],
        "n50": profile["read_length"]["n50"],
        "qv_mean": qv_mean,
        "qv_std": qv_std,
    }, indent=2))


if __name__ == "__main__":
    sys.exit(main())
