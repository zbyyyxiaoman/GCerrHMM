#!/usr/bin/env python3
"""Compare 21-mer overlap between real and simulated FASTQ inputs."""

import argparse
import gzip
import json
import math
import random
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from scipy.stats import pearsonr


def open_fastq(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def reservoir(path, max_reads, k):
    rng = random.Random(20260913)
    sample = []
    for index, record in enumerate(SeqIO.parse(open_fastq(path), "fastq")):
        sequence = str(record.seq).upper()
        if index < max_reads:
            sample.append(sequence)
        else:
            slot = rng.randrange(index + 1)
            if slot < max_reads:
                sample[slot] = sequence

    counts = Counter()
    total = 0
    for sequence in sample:
        for offset in range(len(sequence) - k + 1):
            kmer = sequence[offset:offset + k]
            if "N" not in kmer:
                counts[kmer] += 1
                total += 1
    return counts, total, len(sample)


def compare_counts(real_counts, real_total, sim_counts, sim_total):
    real_set = set(real_counts)
    sim_set = set(sim_counts)
    common = real_set & sim_set
    union = real_set | sim_set
    jaccard = len(common) / len(union) if union else 0.0

    if common:
        real_freq = [real_counts[k] / real_total for k in common]
        sim_freq = [sim_counts[k] / sim_total for k in common]
        correlation = float(pearsonr(real_freq, sim_freq).statistic)
        cosine = float(
            sum(a * b for a, b in zip(real_freq, sim_freq))
            / math.sqrt(sum(a * a for a in real_freq) * sum(b * b for b in sim_freq))
        )
    else:
        correlation = 0.0
        cosine = 0.0

    top_overlap = {}
    for top_n in (100, 1000, 10000):
        real_top = {k for k, _ in real_counts.most_common(top_n)}
        sim_top = {k for k, _ in sim_counts.most_common(top_n)}
        denominator = max(1, min(len(real_top), len(sim_top)))
        top_overlap[str(top_n)] = len(real_top & sim_top) / denominator

    return {
        "real_unique": len(real_set),
        "sim_unique": len(sim_set),
        "common_unique": len(common),
        "jaccard": jaccard,
        "frequency_pearson": correlation,
        "frequency_cosine": cosine,
        "top_overlap_fraction": top_overlap,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--real", required=True)
    parser.add_argument("--sim", required=True)
    parser.add_argument("--k", type=int, default=21)
    parser.add_argument("--reads", type=int, default=250)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    real_counts, real_total, real_reads = reservoir(args.real, args.reads, args.k)
    sim_counts, sim_total, sim_reads = reservoir(args.sim, args.reads, args.k)
    result = {
        "k": args.k,
        "reads_requested": args.reads,
        "reads_analyzed": {"real": real_reads, "simulated": sim_reads},
        "total_kmers": {"real": real_total, "simulated": sim_total},
        **compare_counts(real_counts, real_total, sim_counts, sim_total),
    }
    Path(args.output).write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
