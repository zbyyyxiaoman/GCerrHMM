#!/usr/bin/env python3
"""Measure the error spectrum of a read set against a reference.

Reports global and GC-stratified substitution / insertion / deletion rates
using only primary alignments, so a real read set and a simulated read set can
be compared with the same code. This is the measurement used to choose the
training MAPQ threshold and to re-draw the GC-error curve.

Usage:
  python3 error_spectrum.py --bam aligned.bam --ref ref.fa \
      --label real_ont --output spectrum.json [--window 100] [--bins 10]
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pysam


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, list[str]] = {}
    name = None
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.strip())
    return {key: "".join(value) for key, value in seqs.items()}


def gc_track(sequence: str, window: int) -> np.ndarray:
    """Local GC fraction (0-1) in a centred window, per reference position."""
    is_gc = np.frombuffer(sequence.upper().encode(), dtype=np.uint8)
    is_gc = ((is_gc == ord("G")) | (is_gc == ord("C"))).astype(np.float32)
    kernel = np.ones(window, dtype=np.float32)
    cumulative = np.concatenate([[0.0], np.cumsum(is_gc)])
    half = window // 2
    n = len(sequence)
    left = np.clip(np.arange(n) - half, 0, n)
    right = np.clip(np.arange(n) + half + 1, 0, n)
    counts = cumulative[right] - cumulative[left]
    spans = (right - left).astype(np.float32)
    spans[spans == 0] = 1.0
    return counts / spans


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--label", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--window", type=int, default=100)
    parser.add_argument("--bins", type=int, default=10)
    parser.add_argument("--min-mapq", type=int, default=0)
    args = parser.parse_args()

    reference = read_fasta(Path(args.ref))
    tracks = {name: gc_track(seq, args.window) for name, seq in reference.items()}

    edges = np.linspace(0.0, 1.0, args.bins + 1)
    aligned = np.zeros(args.bins, dtype=np.int64)
    subs = np.zeros(args.bins, dtype=np.int64)
    ins = np.zeros(args.bins, dtype=np.int64)
    dele = np.zeros(args.bins, dtype=np.int64)

    bam = pysam.AlignmentFile(args.bam, "rb")
    reads = 0
    for read in bam:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < args.min_mapq or read.cigar is None:
            continue
        ref_name = read.reference_name
        if ref_name not in reference:
            continue
        reads += 1
        ref_seq = reference[ref_name]
        track = tracks[ref_name]
        query = read.query_sequence
        if query is None:
            continue

        # SAM stores SEQ in the forward-reference orientation, so the stored
        # query_sequence is compared to the forward reference directly. An
        # extra reverse-complement for FLAG 0x10 reads inflates mismatch rates
        # to ~40 % on real data (verified against the same reads without it).
        ref_pos = read.reference_start
        query_pos = 0
        for op, length in read.cigartuples:
            if op in (0, 7, 8):  # M / = / X
                for offset in range(length):
                    rp = ref_pos + offset
                    if rp >= len(ref_seq):
                        break
                    gc = track[rp]
                    b = min(int(gc * args.bins), args.bins - 1)
                    aligned[b] += 1
                    if query_pos + offset < len(query):
                        if query[query_pos + offset].upper() != ref_seq[rp].upper():
                            subs[b] += 1
                ref_pos += length
                query_pos += length
            elif op == 1:  # insertion
                gc = track[min(ref_pos, len(track) - 1)]
                b = min(int(gc * args.bins), args.bins - 1)
                ins[b] += length
                query_pos += length
            elif op == 2:  # deletion
                for offset in range(length):
                    rp = ref_pos + offset
                    if rp >= len(ref_seq):
                        break
                    gc = track[rp]
                    b = min(int(gc * args.bins), args.bins - 1)
                    dele[b] += 1
                ref_pos += length
            elif op == 4:  # soft clip
                query_pos += length
            elif op == 5:  # hard clip
                pass
    bam.close()

    total_aligned = int(aligned.sum())
    total_error = int((subs + ins + dele).sum())
    payload = {
        "label": args.label,
        "bam": str(args.bam),
        "reads": reads,
        "window": args.window,
        "bins": args.bins,
        "global": {
            "aligned_bases": total_aligned,
            "substitution_rate": float(subs.sum() / total_aligned) if total_aligned else 0.0,
            "insertion_rate": float(ins.sum() / total_aligned) if total_aligned else 0.0,
            "deletion_rate": float(dele.sum() / total_aligned) if total_aligned else 0.0,
            "error_rate": float(total_error / total_aligned) if total_aligned else 0.0,
        },
        "by_gc": [],
    }
    for index in range(args.bins):
        denom = int(aligned[index])
        payload["by_gc"].append({
            "gc_low": float(edges[index]),
            "gc_high": float(edges[index + 1]),
            "aligned_bases": denom,
            "substitution_rate": float(subs[index] / denom) if denom else None,
            "insertion_rate": float(ins[index] / denom) if denom else None,
            "deletion_rate": float(dele[index] / denom) if denom else None,
            "error_rate": float((subs[index] + ins[index] + dele[index]) / denom) if denom else None,
        })

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload["global"], indent=2))


if __name__ == "__main__":
    main()
