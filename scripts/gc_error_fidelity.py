#!/usr/bin/env python3
"""Quantify how well a simulated read set reproduces the GC-conditioned error rate.

This is the primary instrument for the GC-conditioning claim. The composite
score contains no term that measures the error rate *conditional on local GC*,
so a curve alone (or a composite) cannot carry the claim.

For one model setting the script reports:

  * Pearson r and Spearman rho between simulated and real per-GC-bin error rates
  * mean absolute deviation (MAD) of the per-bin error rates
  * the per-bin table itself, with support counts

Windows, GC binning and event definition follow the trainer
(`src/train_errhmm.py`): non-overlapping `window_size` bp windows, equal-width
GC bins over the GC fraction, mismatch + insertion base + deletion base over
aligned M bases.

usage:
  python3 gc_error_fidelity.py --ref ref.fa --real-bam real.bam \
      --sim-fastq sim.fastq --sim-bam sim.bam --bins 10 --output out.json
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import numpy as np
import pysam
from scipy import stats


def longest_contig(ref: Path) -> str:
    """Return the name of the longest contig in a FASTA file."""
    import pysam

    with pysam.FastaFile(str(ref)) as fasta:
        return max(
            zip(fasta.references, fasta.lengths),
            key=lambda item: item[1],
        )[0]


def build_pos2bin(ref: pysam.FastaFile, contig: str, window: int,
                  bins: int) -> np.ndarray:
    """Return, for every reference position, the GC bin of its window."""
    seq = ref.fetch(contig).upper().encode()
    arr = np.frombuffer(seq, dtype=np.uint8)
    is_gc = ((arr == ord("G")) | (arr == ord("C"))).astype(np.float64)
    length = arr.size
    n_win = (length + window - 1) // window
    pad = n_win * window - length
    if pad:
        is_gc = np.concatenate([is_gc, np.zeros(pad)])
    frac = is_gc.reshape(n_win, window).mean(axis=1)
    win_bin = np.minimum((frac * bins).astype(np.int32), bins - 1)
    return np.repeat(win_bin.astype(np.int8), window)[:length], seq


def accumulate(bam_path: Path, ref_path: Path, contig: str, window: int,
               bins: int, min_mapq: int = 0) -> dict:
    ref = pysam.FastaFile(str(ref_path))
    pos2bin, ref_seq = build_pos2bin(ref, contig, window, bins)
    ref_arr = np.frombuffer(ref_seq, dtype=np.uint8)
    acgt = np.zeros(256, dtype=bool)
    for ch in b"ACGT":
        acgt[ch] = True

    aligned = np.zeros(bins, dtype=np.int64)
    mismatch = np.zeros(bins, dtype=np.int64)
    ins = np.zeros(bins, dtype=np.int64)
    dele = np.zeros(bins, dtype=np.int64)

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapq:
                continue
            if read.reference_name != contig or read.cigartuples is None:
                continue
            seq = read.query_sequence
            if not seq:
                continue
            read_arr = np.frombuffer(seq.upper().encode(), dtype=np.uint8)
            ref_pos = read.reference_start
            q_pos = 0
            for op, length in read.cigartuples:
                if op in (0, 7, 8):  # M, =, X
                    block_bins = pos2bin[ref_pos:ref_pos + length].astype(np.int64)
                    aligned += np.bincount(block_bins, minlength=bins)
                    rb = ref_arr[ref_pos:ref_pos + length]
                    qb = read_arr[q_pos:q_pos + length]
                    valid = acgt[rb] & acgt[qb]
                    mm = (rb != qb) & valid
                    if mm.any():
                        mismatch += np.bincount(block_bins[mm], minlength=bins)
                    ref_pos += length
                    q_pos += length
                elif op == 1:  # insertion
                    if ref_pos < pos2bin.size:
                        ins[pos2bin[ref_pos]] += length
                    q_pos += length
                elif op == 2:  # deletion
                    b = pos2bin[ref_pos] if ref_pos < pos2bin.size else bins - 1
                    dele[b] += length
                    ref_pos += length
                elif op == 4:  # soft clip
                    q_pos += length
                elif op == 3:  # skipped
                    ref_pos += length
    ref.close()
    return {"aligned": aligned, "mismatch": mismatch, "ins": ins, "del": dele}


def summarise(real: dict, sim: dict, min_support: int = 20000) -> dict:
    with np.errstate(divide="ignore", invalid="ignore"):
        real_rate = np.where(real["aligned"] > 0,
                             (real["mismatch"] + real["ins"] + real["del"])
                             / np.maximum(real["aligned"], 1), np.nan)
        sim_rate = np.where(sim["aligned"] > 0,
                            (sim["mismatch"] + sim["ins"] + sim["del"])
                            / np.maximum(sim["aligned"], 1), np.nan)
    support = (real["aligned"] >= min_support) & (sim["aligned"] >= min_support)
    support &= np.isfinite(real_rate) & np.isfinite(sim_rate)
    used = int(support.sum())
    if used < 3:
        return {"bins_used": used, "pearson_r": None, "spearman_rho": None,
                "mad": None, "per_bin": []}
    real_v = real_rate[support]
    sim_v = sim_rate[support]
    # A constant simulated profile has no defined Pearson correlation. A
    # finite-sample 1-bin run can still vary slightly because its windows have
    # different coverage; such a correlation is not evidence of GC structure.
    if real_v.std() == 0 or sim_v.std() == 0:
        r = float("nan")
        rho = float("nan")
    else:
        r, _ = stats.pearsonr(real_v, sim_v)
        rho, _ = stats.spearmanr(real_v, sim_v)
    mad = float(np.mean(np.abs(real_rate[support] - sim_rate[support])))
    per_bin = []
    for b in range(len(real_rate)):
        per_bin.append({
            "gc_bin": b,
            "real_aligned": int(real["aligned"][b]),
            "sim_aligned": int(sim["aligned"][b]),
            "real_error_rate": (None if not np.isfinite(real_rate[b])
                                else float(real_rate[b])),
            "sim_error_rate": (None if not np.isfinite(sim_rate[b])
                               else float(sim_rate[b])),
            "used": bool(support[b]),
        })
    return {
        "bins_used": used,
        "pearson_r": (None if not np.isfinite(r) else float(r)),
        "spearman_rho": (None if not np.isfinite(rho) else float(rho)),
        "mad": mad,
        "real_mean_error_rate": float(np.mean(real_rate[support])),
        "sim_mean_error_rate": float(np.mean(sim_rate[support])),
        "per_bin": per_bin,
    }


def accumulate_by_window(bam_path: Path, ref_path: Path, contig: str,
                         window: int, min_mapq: int = 0) -> dict:
    """Same events as `accumulate`, but resolved per 100 bp window.

    Aggregating into ten GC bins leaves six to eight usable points, far too few
    for a correlation to mean anything. The window-level arrays keep every
    window that carries aligned bases, so the same comparison runs on thousands
    of points instead of a handful.
    """
    ref = pysam.FastaFile(str(ref_path))
    seq = ref.fetch(contig).upper().encode()
    ref.close()
    ref_arr = np.frombuffer(seq, dtype=np.uint8)
    n_win = (ref_arr.size + window - 1) // window
    aligned = np.zeros(n_win, dtype=np.int64)
    errors = np.zeros(n_win, dtype=np.int64)
    acgt = np.zeros(256, dtype=bool)
    for ch in b"ACGT":
        acgt[ch] = True

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            if read.mapping_quality < min_mapq:
                continue
            if read.reference_name != contig or read.cigartuples is None:
                continue
            seq_r = read.query_sequence
            if not seq_r:
                continue
            read_arr = np.frombuffer(seq_r.upper().encode(), dtype=np.uint8)
            ref_pos, q_pos = read.reference_start, 0
            for op, length in read.cigartuples:
                if op in (0, 7, 8):
                    wins = np.arange(ref_pos, ref_pos + length) // window
                    aligned += np.bincount(wins, minlength=n_win)[:n_win]
                    rb = ref_arr[ref_pos:ref_pos + length]
                    qb = read_arr[q_pos:q_pos + length]
                    mm = (rb != qb) & acgt[rb] & acgt[qb]
                    if mm.any():
                        errors += np.bincount(wins[mm], minlength=n_win)[:n_win]
                    ref_pos += length
                    q_pos += length
                elif op == 1:
                    if ref_pos < ref_arr.size:
                        errors[ref_pos // window] += length
                    q_pos += length
                elif op == 2:
                    if ref_pos < ref_arr.size:
                        errors[ref_pos // window] += length
                    ref_pos += length
                elif op == 4:
                    q_pos += length
                elif op == 3:
                    ref_pos += length
    return {"aligned": aligned, "errors": errors}


def summarise_windows(real: dict, sim: dict, min_support: int = 10,
                      bootstrap: int = 200, seed: int = 11) -> dict:
    """Window-level fidelity plus a bootstrap interval on the correlation."""
    mask = (real["aligned"] >= min_support) & (sim["aligned"] >= min_support)
    n = int(mask.sum())
    if n < 50:
        return {"windows_used": n, "pearson_r": None, "spearman_rho": None,
                "ci_low": None, "ci_high": None, "mad": None}
    real_rate = real["errors"][mask] / real["aligned"][mask]
    sim_rate = sim["errors"][mask] / sim["aligned"][mask]
    r = float(stats.pearsonr(real_rate, sim_rate)[0])
    rho = float(stats.spearmanr(real_rate, sim_rate)[0])
    mad = float(np.mean(np.abs(real_rate - sim_rate)))

    rng = np.random.default_rng(seed)
    idx = np.arange(n)
    draws = []
    for _ in range(bootstrap):
        pick = rng.choice(idx, size=n, replace=True)
        x, y = real_rate[pick], sim_rate[pick]
        if x.std() == 0 or y.std() == 0:
            continue
        draws.append(float(stats.pearsonr(x, y)[0]))
    draws = np.array(draws)
    return {
        "windows_used": n,
        "pearson_r": r,
        "spearman_rho": rho,
        "ci_low": float(np.percentile(draws, 2.5)) if draws.size else None,
        "ci_high": float(np.percentile(draws, 97.5)) if draws.size else None,
        "mad": mad,
        "bootstrap_n": int(draws.size),
        "min_support_bases": min_support,
    }


def window_gc(ref_path: Path, contig: str, window: int) -> np.ndarray:
    """GC fraction of every non-overlapping window of the primary contig."""
    with pysam.FastaFile(str(ref_path)) as fa:
        seq = fa.fetch(contig).upper().encode()
    arr = np.frombuffer(seq, dtype=np.uint8)
    n_win = (arr.size + window - 1) // window
    pad = n_win * window - arr.size
    is_gc = ((arr == ord("G")) | (arr == ord("C"))).astype(np.float64)
    if pad:
        is_gc = np.concatenate([is_gc, np.zeros(pad)])
    return is_gc.reshape(n_win, window).mean(axis=1)


def summarise_curve(real: dict, sim: dict, gc_frac: np.ndarray,
                    curve_bins: int = 50, min_bases: int = 2000,
                    bootstrap: int = 200, seed: int = 11) -> dict:
    """Correlate the real and simulated error curves on ~50 GC points.

    Raw 100 bp windows are far too noisy at 10x (a window holds tens of aligned
    bases), and correlating them attenuates every curve towards zero. Summing
    the windows inside equal-width GC cells keeps thousands of windows in play
    while lifting each point well above the sampling noise, which is what gives
    the correlation real power.
    """
    n_win = min(len(gc_frac), real["aligned"].size, sim["aligned"].size)
    gc_frac = gc_frac[:n_win]
    cell = np.minimum((gc_frac * curve_bins).astype(np.int32), curve_bins - 1)

    def curves(real_al, real_err, sim_al, sim_err):
        r_al = np.bincount(cell, weights=real_al.astype(np.float64),
                           minlength=curve_bins)
        r_er = np.bincount(cell, weights=real_err.astype(np.float64),
                           minlength=curve_bins)
        s_al = np.bincount(cell, weights=sim_al.astype(np.float64),
                           minlength=curve_bins)
        s_er = np.bincount(cell, weights=sim_err.astype(np.float64),
                           minlength=curve_bins)
        ok = (r_al >= min_bases) & (s_al >= min_bases)
        return r_er[ok] / r_al[ok], s_er[ok] / s_al[ok], np.nonzero(ok)[0]

    r_rates, s_rates, cells = curves(real["aligned"], real["errors"],
                                     sim["aligned"], sim["errors"])
    if r_rates.size < 5:
        return {"curve_bins": curve_bins, "points_used": int(r_rates.size),
                "pearson_r": None, "spearman_rho": None, "ci_low": None,
                "ci_high": None, "mad": None, "per_point": []}
    r = float(stats.pearsonr(r_rates, s_rates)[0])
    rho = float(stats.spearmanr(r_rates, s_rates)[0])
    mad = float(np.mean(np.abs(r_rates - s_rates)))

    rng = np.random.default_rng(seed)
    idx = np.arange(n_win)
    draws = []
    for _ in range(bootstrap):
        pick = rng.choice(idx, size=n_win, replace=True)
        rr, ss, _ = curves(real["aligned"][pick], real["errors"][pick],
                           sim["aligned"][pick], sim["errors"][pick])
        if rr.size < 5 or rr.std() == 0 or ss.std() == 0:
            continue
        draws.append(float(stats.pearsonr(rr, ss)[0]))
    draws = np.array(draws)
    return {
        "curve_bins": curve_bins,
        "points_used": int(r_rates.size),
        "min_bases_per_point": min_bases,
        "pearson_r": r,
        "spearman_rho": rho,
        "ci_low": float(np.percentile(draws, 2.5)) if draws.size else None,
        "ci_high": float(np.percentile(draws, 97.5)) if draws.size else None,
        "mad": mad,
        "bootstrap_n": int(draws.size),
        "per_point": [
            {"cell": int(c), "gc_frac": float((c + 0.5) / curve_bins),
             "real_error_rate": float(rr), "sim_error_rate": float(ss)}
            for c, rr, ss in zip(cells, r_rates, s_rates)
        ],
        "real_range": float(r_rates.max() - r_rates.min()),
        "sim_range": float(s_rates.max() - s_rates.min()),
        "real_rate_correlation_with_gc": (
            float(stats.pearsonr(
                np.array([(c + 0.5) / curve_bins for c in cells]), r_rates)[0])
        ),
    }


def ensure_bam(fastq: Path, ref: Path, bam: Path, threads: int) -> None:
    if bam.exists() and bam.stat().st_size > 0:
        return
    bam.parent.mkdir(parents=True, exist_ok=True)
    cmd = (f"minimap2 -ax map-ont -t {threads} {ref} {fastq} 2>/dev/null "
           f"| samtools sort -@ 4 -m 1G -o {bam} -")
    subprocess.run(["bash", "-o", "pipefail", "-c", cmd], check=True)
    subprocess.run(["samtools", "index", str(bam)], check=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", required=True)
    parser.add_argument("--real-bam", required=True)
    parser.add_argument("--sim-bam", required=True)
    parser.add_argument("--sim-fastq", default=None,
                        help="align this FASTQ to --sim-bam when missing")
    parser.add_argument("--contig", default=None,
                        help="primary contig (default: longest in the reference)")
    parser.add_argument("--bins", type=int, default=10)
    parser.add_argument("--window-size", type=int, default=100)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--min-support", type=int, default=20000)
    parser.add_argument("--window-min-support", type=int, default=10,
                        help="aligned bases required in a window, both sides")
    parser.add_argument("--window-bootstrap", type=int, default=200)
    parser.add_argument("--curve-bins", type=int, default=50,
                        help="GC cells used by the higher-power curve analysis")
    parser.add_argument("--curve-min-bases", type=int, default=2000)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    ref_path = Path(args.ref)
    sim_bam = Path(args.sim_bam)
    if not sim_bam.exists() and args.sim_fastq:
        ensure_bam(Path(args.sim_fastq), ref_path, sim_bam, args.threads)
    if not sim_bam.exists():
        raise SystemExit(f"missing simulated BAM: {sim_bam}")

    ref = pysam.FastaFile(str(ref_path))
    contig = args.contig
    if contig is None:
        contig = max(zip(ref.references, ref.lengths), key=lambda x: x[1])[0]
    ref.close()

    real = accumulate(Path(args.real_bam), ref_path, contig, args.window_size,
                      args.bins)
    sim = accumulate(sim_bam, ref_path, contig, args.window_size, args.bins)
    summary = summarise(real, sim, args.min_support)

    # Window-level companion analysis: a correlation over six to eight GC bins
    # carries almost no statistical power, so the same comparison is repeated
    # on every window that has aligned bases in both read sets.
    real_win = accumulate_by_window(Path(args.real_bam), ref_path, contig,
                                    args.window_size)
    sim_win = accumulate_by_window(sim_bam, ref_path, contig,
                                   args.window_size)
    window_level = summarise_windows(real_win, sim_win,
                                     min_support=args.window_min_support,
                                     bootstrap=args.window_bootstrap)
    gc_frac = window_gc(ref_path, contig, args.window_size)
    curve = summarise_curve(real_win, sim_win, gc_frac,
                            curve_bins=args.curve_bins,
                            min_bases=args.curve_min_bases,
                            bootstrap=args.window_bootstrap)
    payload = {
        "reference": str(ref_path),
        "contig": contig,
        "gc_bins": args.bins,
        "window_size": args.window_size,
        "real_bam": str(args.real_bam),
        "sim_bam": str(sim_bam),
        "window_level": window_level,
        "curve": curve,
        **summary,
    }
    Path(args.output).write_text(json.dumps(payload, indent=2))
    print(json.dumps({k: v for k, v in payload.items() if k != "per_bin"},
                     indent=2))


if __name__ == "__main__":
    main()
