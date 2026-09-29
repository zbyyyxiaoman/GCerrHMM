#!/usr/bin/env python3
"""Post-analysis statistics for GC-effect and homopolymer innovation figures."""

import argparse
import csv
import gzip
import itertools
import json
import random
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pysam
from Bio import SeqIO


CANONICAL = set("ACGT")
TOOLS = [
    ("errhmm", "errHMM"),
    ("nanosim", "NanoSim"),
    ("pbsim", "PBSim3"),
    ("badread", "badread"),
]


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def load_reference(path):
    return {
        record.id: str(record.seq).upper()
        for record in SeqIO.parse(str(path), "fasta")
    }


def homopolymer_length(sequence, position):
    if position < 0 or position >= len(sequence):
        return 1
    base = sequence[position]
    if base not in CANONICAL:
        return 1
    left = position
    right = position + 1
    while left > 0 and sequence[left - 1] == base:
        left -= 1
    while right < len(sequence) and sequence[right] == base:
        right += 1
    return right - left


def window_record(gc_windows, contig, ref_pos, window_size):
    start = (ref_pos // window_size) * window_size
    return gc_windows[(contig, start)]


def analyze_bam(bam_path, references, sample_reads, window_size):
    gc_windows = defaultdict(lambda: {
        "errors": 0,
        "aligned_bases": 0,
        "gc_bases": 0,
        "valid_bases": 0,
    })
    homopolymers = defaultdict(lambda: {"deletions": 0, "opportunities": 0})
    reads_seen = 0
    mapped_reads = 0

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for alignment in bam.fetch(until_eof=True):
            if alignment.is_unmapped or alignment.is_secondary or alignment.is_supplementary:
                continue
            if alignment.mapping_quality < 20:
                continue
            query = alignment.query_sequence
            if not query:
                continue
            reads_seen += 1
            if reads_seen > sample_reads:
                break
            mapped_reads += 1
            contig = alignment.reference_name
            reference = references.get(contig)
            if reference is None:
                continue
            query_start = alignment.query_alignment_start
            query_end = alignment.query_alignment_end
            last_ref = None
            for query_pos, ref_pos in alignment.get_aligned_pairs(matches_only=False):
                if query_pos is not None and not (query_start <= query_pos < query_end):
                    # get_aligned_pairs also exposes soft-clipped query bases;
                    # exclude them from the alignment-error denominator.
                    continue
                if ref_pos is not None:
                    last_ref = ref_pos
                    hp_len = min(homopolymer_length(reference, ref_pos), 10)
                    homopolymers[hp_len]["opportunities"] += 1

                if query_pos is not None and ref_pos is not None:
                    query_base = query[query_pos].upper()
                    ref_base = reference[ref_pos]
                    record = window_record(gc_windows, contig, ref_pos, window_size)
                    record["aligned_bases"] += 1
                    if ref_base in CANONICAL:
                        record["valid_bases"] += 1
                    if ref_base in "GC":
                        record["gc_bases"] += 1
                    if query_base != ref_base:
                        record["errors"] += 1
                elif query_pos is not None and last_ref is not None:
                    record = window_record(gc_windows, contig, last_ref, window_size)
                    record["errors"] += 1
                    record["aligned_bases"] += 1
                elif ref_pos is not None:
                    record = window_record(gc_windows, contig, ref_pos, window_size)
                    record["errors"] += 1
                    record["aligned_bases"] += 1
                    hp_len = min(homopolymer_length(reference, ref_pos), 10)
                    homopolymers[hp_len]["deletions"] += 1

    gc_bins = defaultdict(lambda: {"errors": 0, "aligned_bases": 0})
    for record in gc_windows.values():
        if record["valid_bases"] == 0:
            continue
        gc_percent = 100.0 * record["gc_bases"] / record["valid_bases"]
        bin_index = min(9, int(gc_percent // 10))
        gc_bins[bin_index]["errors"] += record["errors"]
        gc_bins[bin_index]["aligned_bases"] += record["aligned_bases"]
    return gc_bins, homopolymers, mapped_reads


def reservoir_reads(path, max_reads, rng):
    sample = []
    for index, record in enumerate(SeqIO.parse(open_text(path), "fastq")):
        sequence = str(record.seq).upper()
        if index < max_reads:
            sample.append(sequence)
        else:
            slot = rng.randrange(index + 1)
            if slot < max_reads:
                sample[slot] = sequence
    return sample


def kmer_vector(sequences, k, kmers):
    counts = Counter()
    total = 0
    for sequence in sequences:
        for offset in range(len(sequence) - k + 1):
            kmer = sequence[offset:offset + k]
            if all(base in CANONICAL for base in kmer):
                counts[kmer] += 1
                total += 1
    if total == 0:
        return np.zeros(len(kmers), dtype=float)
    return np.asarray([counts[kmer] / total for kmer in kmers], dtype=float)


def compute_kmer_pca(datasets, output_csv, output_meta, k, reads_per_sample, samples, seed):
    rng = random.Random(seed)
    max_reads = reads_per_sample * samples
    kmers = ["".join(chars) for chars in itertools.product("ACGT", repeat=k)]
    vectors = []
    labels = []
    for label, path in datasets:
        reads = reservoir_reads(path, max_reads, rng)
        if not reads:
            continue
        for sample_index in range(samples):
            selected = rng.sample(reads, min(reads_per_sample, len(reads)))
            vectors.append(kmer_vector(selected, k, kmers))
            labels.append((label, sample_index))
    matrix = np.vstack(vectors)
    centered = matrix - matrix.mean(axis=0, keepdims=True)
    u, singular_values, _ = np.linalg.svd(centered, full_matrices=False)
    scores = u * singular_values
    variance = singular_values ** 2
    explained = variance / variance.sum()
    with output_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["dataset", "sample", "pc1", "pc2"],
        )
        writer.writeheader()
        for (label, sample_index), point in zip(labels, scores):
            writer.writerow({
                "dataset": label,
                "sample": sample_index,
                "pc1": point[0],
                "pc2": point[1],
            })
    output_meta.write_text(json.dumps({
        "k": k,
        "reads_per_sample": reads_per_sample,
        "samples": samples,
        "explained_variance": explained[:2].tolist(),
    }, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", required=True)
    parser.add_argument("--species", default="Ecoli")
    parser.add_argument("--output-dir")
    parser.add_argument("--reads", type=int, default=20000)
    parser.add_argument("--window-size", type=int, default=100)
    parser.add_argument("--kmer-k", type=int, default=5)
    parser.add_argument("--kmer-reads", type=int, default=100)
    parser.add_argument("--kmer-samples", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20260918)
    parser.add_argument("--skip-kmer", action="store_true")
    args = parser.parse_args()

    project_dir = Path(args.project_dir)
    result_dir = (
        Path(args.output_dir)
        if args.output_dir
        else project_dir / "results" / "innovation"
    )
    bam_dir = project_dir / "results" / "innovation" / "bam"
    result_dir.mkdir(parents=True, exist_ok=True)
    references = load_reference(
        project_dir / "data" / "references" / f"{args.species}_ref.fa"
    )

    bam_inputs = [
        ("real", project_dir / "data" / "real_reads_verified" / f"{args.species}_ont_aligned.bam"),
        *[
            (tool, bam_dir / f"{args.species}_{tool}.sorted.bam")
            for tool, _ in TOOLS
        ],
    ]
    gc_rows = []
    hp_rows = []
    summary_rows = []
    for label, bam_path in bam_inputs:
        if not bam_path.exists():
            summary_rows.append({"dataset": label, "status": "missing", "mapped_reads": 0})
            continue
        gc_bins, homopolymers, mapped_reads = analyze_bam(
            bam_path, references, args.reads, args.window_size
        )
        summary_rows.append({
            "dataset": label,
            "status": "ok",
            "mapped_reads": mapped_reads,
        })
        for bin_index in sorted(gc_bins):
            record = gc_bins[bin_index]
            gc_rows.append({
                "dataset": label,
                "gc_bin": f"{bin_index * 10}-{(bin_index + 1) * 10}",
                "errors": record["errors"],
                "aligned_bases": record["aligned_bases"],
                "error_rate": (
                    record["errors"] / record["aligned_bases"]
                    if record["aligned_bases"] else 0.0
                ),
            })
        for hp_len in sorted(homopolymers):
            record = homopolymers[hp_len]
            hp_rows.append({
                "dataset": label,
                "hp_length": hp_len,
                "deletions": record["deletions"],
                "opportunities": record["opportunities"],
                "deletion_rate": (
                    record["deletions"] / record["opportunities"]
                    if record["opportunities"] else 0.0
                ),
            })

    gc_path = result_dir / "gc_error_curve.csv"
    hp_path = result_dir / "homopolymer_deletion.csv"
    summary_path = result_dir / "bam_summary.csv"
    for path, rows, fields in (
        (gc_path, gc_rows, ["dataset", "gc_bin", "errors", "aligned_bases", "error_rate"]),
        (hp_path, hp_rows, ["dataset", "hp_length", "deletions", "opportunities", "deletion_rate"]),
        (summary_path, summary_rows, ["dataset", "status", "mapped_reads"]),
    ):
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    if not args.skip_kmer:
        data_dir = project_dir / "data"
        kmer_datasets = [
            ("Real", data_dir / "real_reads_verified" / "Ecoli_ont.fastq.gz"),
            ("Route A", data_dir / "simulated" / "route_A_sample" / "Ecoli" / "Ecoli_route_A_sample_r1.fastq.gz"),
            ("Route B", data_dir / "simulated" / "route_B_qshmm" / "Ecoli" / "Ecoli_route_B_qshmm_r1.fastq.gz"),
            ("Route C", data_dir / "simulated" / "route_C_errhmm" / "Ecoli" / "Ecoli_route_C_errhmm_r1.fastq.gz"),
            ("errHMM", data_dir / "simulated" / "cross_tools" / "errhmm" / "Ecoli_fair10k" / "Ecoli_errhmm_10x.fastq"),
            ("NanoSim", data_dir / "simulated" / "cross_tools" / "nanosim" / "Ecoli" / "Ecoli_nanosim_10x.fastq"),
            ("PBSim3", data_dir / "simulated" / "cross_tools" / "pbsim" / "Ecoli" / "Ecoli_pbsim_10x_0001.fq.gz"),
            ("badread", data_dir / "simulated" / "cross_tools" / "badread" / "Ecoli" / "Ecoli_badread_10x.fastq"),
        ]
        compute_kmer_pca(
            kmer_datasets,
            result_dir / "kmer_pca.csv",
            result_dir / "kmer_pca_meta.json",
            args.kmer_k,
            args.kmer_reads,
            args.kmer_samples,
            args.seed,
        )
    print(f"INNOVATION_STATS_DONE={result_dir}")


if __name__ == "__main__":
    main()
