#!/usr/bin/env python3
"""Score an assembled FASTA against a reference in the framework JSON format.

The downstream_tasks.py assembler path always starts from reads. This script
covers the case where the assembly already exists (e.g. an hifiasm run that is
driven by its own wrapper) and only needs the common metrics.

Usage:
  python3 evaluate_assembly_fasta.py --assembly asm.fa --ref ref.fa \
      --threads 4 --output assembly_hifiasm_hifi_Hsapiens_chr21_fair_v2.json
"""

from __future__ import annotations

import argparse
import json
import subprocess
import tempfile
from pathlib import Path


def read_stats(path: Path) -> dict:
    lengths: list[int] = []
    current = 0
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
            else:
                current += len(line.strip())
    if current:
        lengths.append(current)
    lengths.sort(reverse=True)
    total = sum(lengths)
    half = total / 2
    running = 0
    n50 = 0
    for length in lengths:
        running += length
        if running >= half:
            n50 = length
            break
    return {
        "num_contigs": len(lengths),
        "total_length": total,
        "max_contig": lengths[0] if lengths else 0,
        "n50": n50,
        "mean_contig": round(total / len(lengths), 1) if lengths else 0.0,
    }


def alignment_identity(assembly: Path, reference: Path, threads: int) -> dict:
    with tempfile.TemporaryDirectory(prefix="asm_eval_") as tmp:
        paf = Path(tmp) / "asm.paf"
        with paf.open("w") as handle:
            subprocess.run(
                ["minimap2", "-cx", "asm5", "-t", str(threads),
                 str(reference), str(assembly)],
                stdout=handle, stderr=subprocess.DEVNULL, check=True,
            )
        aligned = 0
        matched = 0
        for line in paf.read_text().splitlines():
            fields = line.split("\t")
            if len(fields) < 12:
                continue
            aligned += int(fields[10])
            matched += int(fields[9])
    return {
        "total_aligned_bases": aligned,
        "matched_bases": matched,
        "reference_identity": round(matched / aligned, 4) if aligned else 0.0,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    assembly = Path(args.assembly)
    reference = Path(args.ref)
    if not assembly.exists():
        raise SystemExit(f"missing assembly: {assembly}")
    if not reference.exists():
        raise SystemExit(f"missing reference: {reference}")

    payload = {"assembler_source": str(assembly)}
    payload.update(read_stats(assembly))
    payload.update(alignment_identity(assembly, reference, args.threads))

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
