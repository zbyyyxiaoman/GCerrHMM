#!/usr/bin/env python3
"""Compute only the pre-specified 10-bin GC-error curve endpoint.

This is a numerical fast path for the confirmation cohort. It uses the same
100-bp windows, equal-width GC bins and CIGAR event definitions as
``gc_error_fidelity.py``, but skips the descriptive window-level bootstrap.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pysam
from scipy import stats


def build_pos2bin(
    fasta: pysam.FastaFile,
    contig: str,
    window: int,
    bins: int,
) -> tuple[np.ndarray, bytes]:
    sequence = fasta.fetch(contig).upper().encode()
    ref = np.frombuffer(sequence, dtype=np.uint8)
    is_gc = ((ref == ord("G")) | (ref == ord("C"))).astype(np.float64)
    length = ref.size
    n_windows = (length + window - 1) // window
    padding = n_windows * window - length
    if padding:
        is_gc = np.concatenate([is_gc, np.zeros(padding)])
    fraction = is_gc.reshape(n_windows, window).mean(axis=1)
    window_bin = np.minimum((fraction * bins).astype(np.int32), bins - 1)
    return np.repeat(window_bin.astype(np.int8), window)[:length], sequence


def accumulate(
    bam_path: Path,
    ref_path: Path,
    contig: str,
    window: int,
    bins: int,
    min_mapq: int = 0,
) -> dict[str, np.ndarray]:
    with pysam.FastaFile(str(ref_path)) as fasta:
        pos2bin, ref_sequence = build_pos2bin(fasta, contig, window, bins)
    ref = np.frombuffer(ref_sequence, dtype=np.uint8)
    acgt = np.zeros(256, dtype=bool)
    for base in b"ACGT":
        acgt[base] = True

    aligned = np.zeros(bins, dtype=np.int64)
    mismatch = np.zeros(bins, dtype=np.int64)
    insertions = np.zeros(bins, dtype=np.int64)
    deletions = np.zeros(bins, dtype=np.int64)

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapq:
                continue
            if read.reference_name != contig or read.cigartuples is None:
                continue
            sequence = read.query_sequence
            if not sequence:
                continue
            read_bases = np.frombuffer(sequence.upper().encode(), dtype=np.uint8)
            ref_pos = read.reference_start
            query_pos = 0
            for operation, length in read.cigartuples:
                if operation in (0, 7, 8):
                    block_bins = pos2bin[ref_pos:ref_pos + length].astype(np.int64)
                    aligned += np.bincount(block_bins, minlength=bins)
                    ref_block = ref[ref_pos:ref_pos + length]
                    query_block = read_bases[query_pos:query_pos + length]
                    valid = acgt[ref_block] & acgt[query_block]
                    substitutions = (ref_block != query_block) & valid
                    if substitutions.any():
                        mismatch += np.bincount(
                            block_bins[substitutions],
                            minlength=bins,
                        )
                    ref_pos += length
                    query_pos += length
                elif operation == 1:
                    if ref_pos < pos2bin.size:
                        insertions[pos2bin[ref_pos]] += length
                    query_pos += length
                elif operation == 2:
                    index = pos2bin[ref_pos] if ref_pos < pos2bin.size else bins - 1
                    deletions[index] += length
                    ref_pos += length
                elif operation == 4:
                    query_pos += length
                elif operation == 3:
                    ref_pos += length
    return {
        "aligned": aligned,
        "mismatch": mismatch,
        "insertions": insertions,
        "deletions": deletions,
    }


def summarise(
    real: dict[str, np.ndarray],
    simulated: dict[str, np.ndarray],
    min_support: int,
) -> dict:
    with np.errstate(divide="ignore", invalid="ignore"):
        real_rate = np.where(
            real["aligned"] > 0,
            (
                real["mismatch"] + real["insertions"] + real["deletions"]
            )
            / np.maximum(real["aligned"], 1),
            np.nan,
        )
        simulated_rate = np.where(
            simulated["aligned"] > 0,
            (
                simulated["mismatch"]
                + simulated["insertions"]
                + simulated["deletions"]
            )
            / np.maximum(simulated["aligned"], 1),
            np.nan,
        )
    support = (
        (real["aligned"] >= min_support)
        & (simulated["aligned"] >= min_support)
        & np.isfinite(real_rate)
        & np.isfinite(simulated_rate)
    )
    used = int(support.sum())
    if used < 3:
        return {
            "bins_used": used,
            "pearson_r": None,
            "spearman_rho": None,
            "mad": None,
        }
    real_values = real_rate[support]
    simulated_values = simulated_rate[support]
    if real_values.std() == 0 or simulated_values.std() == 0:
        pearson = float("nan")
        spearman = float("nan")
    else:
        pearson = float(stats.pearsonr(real_values, simulated_values)[0])
        spearman = float(stats.spearmanr(real_values, simulated_values)[0])
    return {
        "bins_used": used,
        "pearson_r": None if not np.isfinite(pearson) else pearson,
        "spearman_rho": None if not np.isfinite(spearman) else spearman,
        "mad": float(np.mean(np.abs(real_values - simulated_values))),
        "real_mean_error_rate": float(np.mean(real_values)),
        "sim_mean_error_rate": float(np.mean(simulated_values)),
        "real_counts": {
            key: value.astype(int).tolist() for key, value in real.items()
        },
        "simulated_counts": {
            key: value.astype(int).tolist() for key, value in simulated.items()
        },
    }


def save_counts(path: Path, counts: dict[str, np.ndarray]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez(path, **counts)


def load_counts(path: Path) -> dict[str, np.ndarray]:
    with np.load(path) as payload:
        return {key: payload[key] for key in payload.files}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", type=Path, required=True)
    parser.add_argument("--real-bam", type=Path, required=True)
    parser.add_argument("--sim-bam", type=Path, required=True)
    parser.add_argument(
        "--real-counts",
        type=Path,
        help="Reuse or write cached real-read counts at this path.",
    )
    parser.add_argument("--contig")
    parser.add_argument("--bins", type=int, default=10)
    parser.add_argument("--window-size", type=int, default=100)
    parser.add_argument("--min-support", type=int, default=20000)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with pysam.FastaFile(str(args.ref)) as fasta:
        contig = args.contig or max(
            zip(fasta.references, fasta.lengths),
            key=lambda item: item[1],
        )[0]

    if args.real_counts and args.real_counts.exists():
        real = load_counts(args.real_counts)
    else:
        real = accumulate(
            args.real_bam,
            args.ref,
            contig,
            args.window_size,
            args.bins,
        )
        if args.real_counts:
            save_counts(args.real_counts, real)
    simulated = accumulate(
        args.sim_bam,
        args.ref,
        contig,
        args.window_size,
        args.bins,
    )
    payload = {
        "reference": str(args.ref),
        "real_bam": str(args.real_bam),
        "sim_bam": str(args.sim_bam),
        "contig": contig,
        "gc_bins": args.bins,
        "window_size": args.window_size,
        "endpoint_definition": "10-bin GC-stratified curve Pearson r",
        "real_counts_cache": (
            None if args.real_counts is None else str(args.real_counts)
        ),
        **summarise(real, simulated, args.min_support),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({key: value for key, value in payload.items()
                      if not key.endswith("_counts")}, indent=2))


if __name__ == "__main__":
    main()
