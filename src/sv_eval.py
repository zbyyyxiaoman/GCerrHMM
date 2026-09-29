#!/usr/bin/env python3
"""Structural-variant calling and truth comparison for simulated reads.

Maps simulated long reads, calls SVs with Sniffles2, and compares the call
set against the spike-in truth VCF produced by ``sv_spikein.py`` using the
same caller and parameters for every simulator.

Usage:
  python3 sv_eval.py --ref ref.fa --truth truth_sv.vcf --sim sim.fastq \
      --threads 4 --output variant_sv_<tool>_<species>_sv.json
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import shutil
import subprocess
import tempfile
from pathlib import Path

MIN_SVLEN = 50
INS_POS_TOL = 300
OVERLAP_MIN = 0.5


def run(cmd: str, check: bool = True) -> subprocess.CompletedProcess:
    print(f"  [CMD] {cmd[:140]}", flush=True)
    res = subprocess.run(
        ["bash", "-o", "pipefail", "-c", cmd], capture_output=True, text=True
    )
    if check and res.returncode != 0:
        raise RuntimeError(f"command failed: {cmd}\n{res.stderr[-2000:]}")
    return res


def open_text(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def parse_info(info: str) -> dict[str, str]:
    out = {}
    for item in info.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        else:
            out[item] = "True"
    return out


def read_sv_vcf(path: Path) -> list[dict]:
    """Return SV records with chrom/pos/end/svtype/svlen."""
    records: list[dict] = []
    if not path.exists():
        return records
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            chrom, pos = f[0], int(f[1])
            info = parse_info(f[7])
            svtype = info.get("SVTYPE", "")
            if not svtype:
                alt = f[4]
                if alt.startswith("<") and alt.endswith(">"):
                    svtype = alt[1:-1]
            svtype = svtype.upper()
            if svtype in {"BND", "TRA", "CNV"}:
                continue
            end = int(info.get("END", pos))
            svlen = info.get("SVLEN", "")
            try:
                svlen_i = int(str(svlen).split(",")[0])
            except (ValueError, IndexError):
                svlen_i = end - pos if svtype != "INS" else 0
            if abs(svlen_i) < MIN_SVLEN:
                continue
            records.append(
                {
                    "chrom": chrom,
                    "pos": pos,
                    "end": end if end >= pos else pos,
                    "svtype": svtype,
                    "svlen": abs(svlen_i),
                }
            )
    return records


def overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> float:
    lo = max(a_start, b_start)
    hi = min(a_end, b_end)
    inter = max(0, hi - lo)
    if inter == 0:
        return 0.0
    a_len = max(1, a_end - a_start)
    b_len = max(1, b_end - b_start)
    return inter / min(a_len, b_len)


def match(truth: list[dict], calls: list[dict]) -> tuple[list, list, list]:
    """Greedy one-to-one matching by best score within each SV type."""
    used_calls: set[int] = set()
    pairs: list[tuple[int, int, float]] = []
    for ti, t in enumerate(truth):
        for ci, c in enumerate(calls):
            if c["chrom"] != t["chrom"] or c["svtype"] != t["svtype"]:
                continue
            if t["svtype"] == "INS":
                if abs(c["pos"] - t["pos"]) > INS_POS_TOL:
                    continue
                len_ok = abs(c["svlen"] - t["svlen"]) <= max(100, 0.5 * t["svlen"])
                if not len_ok:
                    continue
                score = 1.0 - abs(c["pos"] - t["pos"]) / (INS_POS_TOL * 4)
            else:
                ov = overlap(t["pos"], t["end"], c["pos"], c["end"])
                if ov < OVERLAP_MIN:
                    continue
                score = ov
            pairs.append((ti, ci, score))

    pairs.sort(key=lambda p: p[2], reverse=True)
    matched_truth: set[int] = set()
    tp_pairs: list[tuple[int, int]] = []
    for ti, ci, _ in pairs:
        if ti in matched_truth or ci in used_calls:
            continue
        matched_truth.add(ti)
        used_calls.add(ci)
        tp_pairs.append((ti, ci))

    fp = [c for ci, c in enumerate(calls) if ci not in used_calls]
    fn = [t for ti, t in enumerate(truth) if ti not in matched_truth]
    return tp_pairs, fp, fn


def prf(tp: int, fp: int, fn: int) -> dict:
    precision = tp / (tp + fp) if tp + fp else 0.0
    recall = tp / (tp + fn) if tp + fn else 0.0
    f1 = (
        2 * precision * recall / (precision + recall)
        if precision + recall
        else 0.0
    )
    return {
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "precision": round(precision, 4),
        "recall": round(recall, 4),
        "f1": round(f1, 4),
        "precision_ci95": wilson(tp, tp + fp),
        "recall_ci95": wilson(tp, tp + fn),
    }


def wilson(successes: int, total: int, z: float = 1.959963985) -> list[float]:
    """Wilson score interval, which behaves sanely for small event counts."""
    if total <= 0:
        return [0.0, 0.0]
    phat = successes / total
    denom = 1 + z * z / total
    centre = phat + z * z / (2 * total)
    margin = z * math.sqrt(phat * (1 - phat) / total + z * z / (4 * total * total))
    return [
        round(max(0.0, (centre - margin) / denom), 4),
        round(min(1.0, (centre + margin) / denom), 4),
    ]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--truth", required=True, help="SV truth VCF")
    ap.add_argument("--sim", required=True, help="simulated reads (fastq/fastq.gz)")
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--output", required=True)
    ap.add_argument("--workdir", default=None)
    args = ap.parse_args()

    threads = args.threads
    work = Path(args.workdir) if args.workdir else Path(tempfile.mkdtemp(prefix="sv_eval_"))
    work.mkdir(parents=True, exist_ok=True)
    bam = work / "aligned.bam"
    vcf = work / "sv_calls.vcf"

    run(
        f"minimap2 -ax map-ont --secondary=no -t {threads} {args.ref} {args.sim} "
        f"| samtools sort -@ {threads} -o {bam}"
    )
    run(f"samtools index {bam}")
    run(
        f"sniffles -i {bam} -v {vcf} --threads {threads} "
        f"--minsvlen {MIN_SVLEN} --mapq 20"
    )

    truth = read_sv_vcf(Path(args.truth))
    calls = read_sv_vcf(vcf)
    tp_pairs, fp, fn = match(truth, calls)

    types = sorted({t["svtype"] for t in truth} | {c["svtype"] for c in calls})
    per_type = {}
    for svtype in types:
        t_sub = [t for t in truth if t["svtype"] == svtype]
        c_sub = [c for c in calls if c["svtype"] == svtype]
        tp_sub, fp_sub, fn_sub = match(t_sub, c_sub)
        per_type[svtype] = prf(len(tp_sub), len(fp_sub), len(fn_sub))

    size_bins = {"50-500": (50, 500), "500-2000": (500, 2000), "2000+": (2000, 10**9)}
    per_size = {}
    for label, (lo, hi) in size_bins.items():
        t_sub = [t for t in truth if lo <= t["svlen"] < hi]
        c_sub = [c for c in calls if lo <= c["svlen"] < hi]
        tp_sub, fp_sub, fn_sub = match(t_sub, c_sub)
        per_size[label] = prf(len(tp_sub), len(fp_sub), len(fn_sub))

    result = {
        "sv_caller": "sniffles2",
        "min_svlen": MIN_SVLEN,
        "mapq_min": 20,
        "truth_total": len(truth),
        "called_total": len(calls),
        "overall": prf(len(tp_pairs), len(fp), len(fn)),
        "by_type": per_type,
        "by_size": per_size,
        "truth_counts": {t: sum(1 for x in truth if x["svtype"] == t) for t in types},
        "call_counts": {t: sum(1 for x in calls if x["svtype"] == t) for t in types},
    }

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        json.dump(result, fh, indent=2, sort_keys=True)
    print(json.dumps(result, indent=2, sort_keys=True))

    if not args.workdir:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
