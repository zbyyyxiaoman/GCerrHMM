#!/usr/bin/env python3
"""Create a deterministic variant truth set and a variant-bearing reference.

The generated FASTA is used by Phase 3 read simulation so simulated reads
contain real germline-like variants.  The VCF is kept against the original
reference and is passed to downstream_tasks.py variant evaluation.
"""

import argparse
import json
import random
from collections import Counter
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


BASES = set("ACGT")


def _clean_base(base):
    return base.upper() if base.upper() in BASES else None


def _generate_variants(record, rng, snp_rate, indel_rate, min_distance):
    seq = str(record.seq).upper()
    length = len(seq)
    variants = []

    n_snps = int(round(length * snp_rate))
    n_indels = int(round(length * indel_rate))
    if n_snps > 0:
        for pos in rng.sample(range(1, length + 1), min(n_snps, length)):
            ref = _clean_base(seq[pos - 1])
            if ref is None:
                continue
            alt = rng.choice([b for b in "ACGT" if b != ref])
            variants.append({"type": "snp", "pos": pos, "ref": ref, "alt": alt})

    if n_indels > 0:
        positions = rng.sample(range(1, length + 1), min(n_indels, length))
        for pos in positions:
            ref = _clean_base(seq[pos - 1])
            if ref is None:
                continue
            if rng.random() < 0.5:
                ins_len = rng.choice([1, 1, 1, 2, 2, 3])
                ins = "".join(rng.choice("ACGT") for _ in range(ins_len))
                variants.append({
                    "type": "ins",
                    "pos": pos,
                    "ref": ref,
                    "alt": ref + ins,
                    "ins": ins,
                    "length": ins_len,
                })
            else:
                del_len = rng.choice([1, 1, 1, 2, 2, 3])
                end = pos + del_len
                if end > length:
                    continue
                ref_seq = seq[pos - 1:end]
                if any(_clean_base(b) is None for b in ref_seq):
                    continue
                variants.append({
                    "type": "del",
                    "pos": pos,
                    "ref": ref_seq,
                    "alt": ref_seq[0],
                    "length": del_len,
                    "ref_len": len(ref_seq),
                })

    # Keep variants separated enough that they cannot overlap or create
    # ambiguous truth records during left normalization.
    variants.sort(key=lambda v: (v["pos"], v["type"]))
    kept = []
    last_end = -10**9
    for v in variants:
        start = v["pos"]
        ref_len = v.get("ref_len", v.get("length", 1) + 1 if v["type"] == "del" else 1)
        end = start + ref_len - 1
        if start - last_end >= min_distance:
            kept.append(v)
            last_end = end
    return kept


def _apply_variants(seq, variants):
    chunks = []
    cursor = 0
    for v in variants:
        p = v["pos"] - 1
        if p < cursor:
            continue
        if v["type"] == "snp":
            chunks.append(seq[cursor:p])
            chunks.append(v["alt"])
            cursor = p + 1
        elif v["type"] == "ins":
            chunks.append(seq[cursor:p + 1])
            chunks.append(v["ins"])
            cursor = p + 1
        elif v["type"] == "del":
            chunks.append(seq[cursor:p])
            chunks.append(seq[p])
            cursor = p + v.get("ref_len", v["length"] + 1)
    chunks.append(seq[cursor:])
    return "".join(chunks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref", required=True, help="Original reference FASTA")
    parser.add_argument("--output-fasta", required=True, help="Variant-bearing FASTA")
    parser.add_argument("--output-vcf", required=True, help="Truth VCF against original reference")
    parser.add_argument("--output-stats", help="Optional JSON stats path")
    parser.add_argument("--snp-rate", type=float, default=0.001)
    parser.add_argument("--indel-rate", type=float, default=0.0001)
    parser.add_argument("--min-distance", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    records = list(SeqIO.parse(args.ref, "fasta"))
    if not records:
        raise SystemExit(f"No FASTA records found in {args.ref}")

    out_fasta = Path(args.output_fasta)
    out_vcf = Path(args.output_vcf)
    out_fasta.parent.mkdir(parents=True, exist_ok=True)
    out_vcf.parent.mkdir(parents=True, exist_ok=True)

    all_variants = []
    output_records = []
    genome_length = 0

    for record in records:
        seq = str(record.seq).upper()
        genome_length += len(seq)
        variants = _generate_variants(
            record,
            rng,
            args.snp_rate,
            args.indel_rate,
            args.min_distance,
        )
        for v in variants:
            v["chrom"] = record.id
        all_variants.extend(variants)
        mutated = _apply_variants(seq, variants)
        output_records.append(SeqRecord(Seq(mutated), id=record.id, description=record.description))

    with open(out_fasta, "w", encoding="utf-8") as fh:
        SeqIO.write(output_records, fh, "fasta")

    contig_lengths = {r.id: len(str(r.seq)) for r in records}
    all_variants.sort(key=lambda v: (v["chrom"], v["pos"]))

    with open(out_vcf, "w", encoding="utf-8") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write(f"##source=create_truth_variants.py,seed={args.seed},snp_rate={args.snp_rate},indel_rate={args.indel_rate}\n")
        for chrom, length in contig_lengths.items():
            fh.write(f"##contig=<ID={chrom},length={length}>\n")
        fh.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Variant type">\n')
        fh.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of variant">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for i, v in enumerate(all_variants, start=1):
            svtype = v["type"].upper()
            ref_len = v.get("ref_len", v.get("length", 1) + 1 if v["type"] == "del" else 1)
            end = v["pos"] + ref_len - 1
            info = f"SVTYPE={svtype};END={end}"
            fh.write(
                f"{v['chrom']}\t{v['pos']}\ttruth_{i}\t{v['ref']}\t{v['alt']}\t.\tPASS\t{info}\n"
            )

    snp_count = sum(1 for v in all_variants if v["type"] == "snp")
    indel_count = len(all_variants) - snp_count
    ins_count = sum(1 for v in all_variants if v["type"] == "ins")
    del_count = sum(1 for v in all_variants if v["type"] == "del")
    indel_lengths = Counter(v.get("length", 1) for v in all_variants if v["type"] != "snp")
    stats = {
        "reference": str(Path(args.ref)),
        "output_fasta": str(out_fasta),
        "output_vcf": str(out_vcf),
        "genome_length": genome_length,
        "snp_rate": args.snp_rate,
        "indel_rate": args.indel_rate,
        "min_distance": args.min_distance,
        "seed": args.seed,
        "snp_count": snp_count,
        "indel_count": indel_count,
        "insertion_count": ins_count,
        "deletion_count": del_count,
        "indel_length_distribution": dict(sorted(indel_lengths.items())),
        "total_variants": len(all_variants),
    }
    stats_path = Path(args.output_stats or str(out_vcf).replace(".vcf", ".stats.json"))
    with open(stats_path, "w", encoding="utf-8") as fh:
        json.dump(stats, fh, indent=2)

    print(json.dumps(stats, indent=2))


if __name__ == "__main__":
    main()
