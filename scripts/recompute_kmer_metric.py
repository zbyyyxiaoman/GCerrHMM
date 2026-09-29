#!/usr/bin/env python3
"""Recompute the k-mer sub-score with more reads and a bootstrap interval.

The composite's k-mer term correlates the real and simulated k-mer frequency
vectors over the (tiny) intersection of their k-mer sets, using 250 reads. That
intersection holds a few hundred k-mers, so the correlation moves by up to
16.7 points between replicates of the *same* route.

This script keeps the statistic identical - same k, same Pearson correlation
over the intersection of the per-side top-K maps - and only

  1. samples more reads (default 1000 instead of 250), and
  2. resamples reads with replacement to attach a 95% bootstrap interval.

Sampling and interval estimation reduce variance; they do not change the
estimator, so this is not metric shopping.

usage:
  python3 recompute_kmer_metric.py --real real.fastq.gz --sim sim.fastq \
      --reads 1000 --bootstrap 50 --output kmer.json
"""

from __future__ import annotations

import argparse
import gzip
import json
import random
from pathlib import Path

import numpy as np
from scipy import stats


def opener(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def read_sequences(path: Path, limit: int) -> list[str]:
    seqs = []
    with opener(path) as handle:
        for i, line in enumerate(handle):
            if i % 4 == 1:
                seqs.append(line.strip())
                if len(seqs) >= limit:
                    break
    return seqs


CODE = np.full(256, 255, dtype=np.uint8)
for _base, _val in zip(b"ACGT", range(4)):
    CODE[_base] = _val


def kmers(seq: str, k: int) -> np.ndarray:
    a = np.frombuffer(seq.upper().encode(), dtype=np.uint8)
    if a.size < k:
        return np.empty(0, dtype=np.int64)
    c = CODE[a]
    n = c.size - k + 1
    val = np.zeros(n, dtype=np.int64)
    bad = np.zeros(n, dtype=bool)
    for j in range(k):
        chunk = c[j:j + n]
        val = (val << 2) | chunk
        bad |= chunk == 255
    return val[~bad]


def per_read_counts(seqs: list[str], keys: np.ndarray, k: int):
    """Per-read (key index, count) pairs restricted to `keys`."""
    lookup = {int(key): idx for idx, key in enumerate(keys)}
    per_read = []
    total = 0
    for seq in seqs:
        vals = kmers(seq, k)
        if vals.size == 0:
            per_read.append((np.empty(0, dtype=np.int64),
                             np.empty(0, dtype=np.int64)))
            continue
        uniq, counts = np.unique(vals, return_counts=True)
        idx = np.array([lookup.get(int(u), -1) for u in uniq], dtype=np.int64)
        keep = idx >= 0
        per_read.append((idx[keep], counts[keep].astype(np.int64)))
        total += int(counts.sum())
    return per_read, total


def vector_from(per_read, picks, n_keys: int) -> np.ndarray:
    vec = np.zeros(n_keys, dtype=np.float64)
    for i in picks:
        idx, cnt = per_read[i]
        if idx.size:
            np.add.at(vec, idx, cnt)
    return vec


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--real", required=True)
    parser.add_argument("--sim", required=True)
    parser.add_argument("--reads", type=int, default=1000)
    parser.add_argument("--k", type=int, default=21)
    parser.add_argument("--top", type=int, default=10000)
    parser.add_argument("--bootstrap", type=int, default=50)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    real_seqs = read_sequences(Path(args.real), args.reads)
    sim_seqs = read_sequences(Path(args.sim), args.reads)
    n = min(len(real_seqs), len(sim_seqs))
    real_seqs, sim_seqs = real_seqs[:n], sim_seqs[:n]

    def top_keys(seqs):
        vals = np.concatenate([kmers(s, args.k) for s in seqs]) \
            if seqs else np.empty(0, dtype=np.int64)
        uniq, counts = np.unique(vals, return_counts=True)
        order = np.argsort(counts)[::-1][:args.top]
        return uniq[order]

    real_keys = top_keys(real_seqs)
    sim_keys = top_keys(sim_seqs)
    common = np.intersect1d(real_keys, sim_keys)
    if common.size < 10:
        raise SystemExit(f"too few shared k-mers ({common.size})")

    real_reads, _ = per_read_counts(real_seqs, common, args.k)
    sim_reads, _ = per_read_counts(sim_seqs, common, args.k)

    all_picks = np.arange(n)
    real_vec = vector_from(real_reads, all_picks, common.size)
    sim_vec = vector_from(sim_reads, all_picks, common.size)
    r_point, _ = stats.pearsonr(real_vec, sim_vec)

    rng = random.Random(args.seed)
    boot = []
    for _ in range(args.bootstrap):
        picks = [rng.randrange(n) for _ in range(n)]
        rv = vector_from(real_reads, picks, common.size)
        sv = vector_from(sim_reads, picks, common.size)
        if rv.std() == 0 or sv.std() == 0:
            continue
        boot.append(float(stats.pearsonr(rv, sv)[0]))

    payload = {
        "statistic": "pearson_r_over_shared_kmers",
        "k": args.k,
        "reads_sampled": n,
        "shared_kmers": int(common.size),
        "pearson_r": float(r_point),
        "bootstrap_n": len(boot),
        "bootstrap_ci_low": float(np.percentile(boot, 2.5)) if boot else None,
        "bootstrap_ci_high": float(np.percentile(boot, 97.5)) if boot else None,
        "bootstrap_sd": float(np.std(boot)) if boot else None,
        "real": str(args.real),
        "sim": str(args.sim),
    }
    Path(args.output).write_text(json.dumps(payload, indent=2))
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
