#!/usr/bin/env python3
"""Build an SV spike-in donor reference and its truth VCF.

The project's truth sets (GIAB / SGD-derived panels) contain only small
variants, so structural-variant (SV) calling cannot be evaluated on them.
This module injects a deterministic SV panel into the existing variant
reference (the donor genome used for simulation) and records the events in
the coordinate system of the original reference, which is what read mappers
and SV callers report against.

Usage:
  python3 sv_spikein.py \
      --ref data/references/Ecoli_ref.fa \
      --variant-ref data/truth/Ecoli_variant_ref.fa \
      --truth-vcf data/truth/Ecoli_truth.vcf \
      --output-ref data/truth/sv/Ecoli_sv_donor.fa \
      --output-vcf data/truth/sv/Ecoli_sv_truth.vcf
"""

from __future__ import annotations

import argparse
import gzip
import random
import re
from pathlib import Path

EVENT_SIZES = {
    "DEL": (200, 1000, 5000, 20000),
    "INS": (200, 1000, 5000),
    "INV": (500, 2000, 10000),
    "DUP": (500, 1000, 3000),
}

EVENTS_PER_COMBO = 3
MIN_CONTIG = 200_000
MARGIN = 50_000

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def read_fasta(path: Path) -> list[tuple[str, str]]:
    seqs: list[tuple[str, str]] = []
    name = None
    chunks: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(chunks)))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs.append((name, "".join(chunks)))
    return seqs


def write_fasta(path: Path, seqs: list[tuple[str, str]], width: int = 60) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for name, seq in seqs:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")


def open_maybe_gzip(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path)


def parse_variants(path: Path) -> dict[str, list[tuple[int, str, str]]]:
    """Return {contig: [(pos1, ref, alt), ...]} sorted by position."""
    out: dict[str, list[tuple[int, str, str]]] = {}
    pattern = re.compile(r"[ACGTNacgtn]+")
    with open_maybe_gzip(path) as fh:
        for line in fh:
            if not line or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            chrom, pos, ref, alt = f[0], int(f[1]), f[3], f[4]
            if alt.startswith("<"):
                continue
            if not pattern.fullmatch(alt) or not pattern.fullmatch(ref):
                continue
            out.setdefault(chrom, []).append((pos, ref.upper(), alt.upper()))
    for chrom in out:
        out[chrom].sort()
    return out


def build_offset_map(
    variants: dict[str, list[tuple[int, str, str]]],
) -> dict[str, list[tuple[int, int]]]:
    """Cumulative length delta at each variant position.

    Mirrors ``bcftools consensus``: variants are applied left to right and
    each one shifts all downstream coordinates by len(ALT) - len(REF).
    """
    maps: dict[str, list[tuple[int, int]]] = {}
    for chrom, rows in variants.items():
        entries: list[tuple[int, int]] = []
        delta = 0
        for pos, ref, alt in rows:
            delta += len(alt) - len(ref)
            entries.append((pos, delta))
        maps[chrom] = entries
    return maps


def to_variant_coord(entries: list[tuple[int, int]], pos: int) -> int:
    """Convert an original-reference position into variant-reference space."""
    delta = 0
    for vpos, vdelta in entries:
        if vpos <= pos:
            delta = vdelta
        else:
            break
    return pos + delta


def random_sequence(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


def plan_events(
    lengths: dict[str, int], rng: random.Random, events_per_combo: int = EVENTS_PER_COMBO
) -> dict[str, list[dict]]:
    """Spread one event of each (type, size) repeatedly across the genome."""
    combos = [
        (svtype, size)
        for svtype, sizes in EVENT_SIZES.items()
        for size in sizes
        for _ in range(events_per_combo)
    ]
    usable = {n: L for n, L in lengths.items() if L >= MIN_CONTIG} or dict(lengths)
    total = sum(usable.values())

    slots: list[tuple[str, int]] = []
    for name, length in sorted(usable.items()):
        n = max(1, round(len(combos) * length / total))
        span = max(1, length - 2 * MARGIN)
        for i in range(n):
            pos = MARGIN + int((i + 0.5) * span / n) + rng.randint(-5_000, 5_000)
            pos = max(1, min(length - 1, pos))
            slots.append((name, pos))

    rng.shuffle(combos)
    events: dict[str, list[dict]] = {}
    for (svtype, size), (name, pos) in zip(combos, slots):
        events.setdefault(name, []).append(
            {"type": svtype, "orig_pos": pos, "pos": pos, "size": size}
        )
    return events


def inject(seq: str, events: list[dict], rng: random.Random) -> tuple[str, list[dict]]:
    """Apply events to one contig and emit truth records.

    ``events`` carry both the injection coordinate (``pos``) and the truth
    coordinate in the original reference (``orig_pos``). Events are applied
    from the highest coordinate downwards so earlier coordinates stay valid.
    """
    truth: list[dict] = []
    for ev in sorted(events, key=lambda e: e["pos"], reverse=True):
        pos = ev["pos"]
        size = ev["size"]
        svtype = ev["type"]
        start = pos - 1
        end = start + size
        if start < 0 or end >= len(seq):
            continue
        orig = ev["orig_pos"]
        if svtype == "DEL":
            seq = seq[:start] + seq[end:]
            truth.append(
                {"pos": orig, "end": orig + size, "svtype": "DEL", "svlen": -size}
            )
        elif svtype == "INS":
            seq = seq[:start] + random_sequence(rng, size) + seq[start:]
            truth.append({"pos": orig, "end": orig, "svtype": "INS", "svlen": size})
        elif svtype == "INV":
            seq = seq[:start] + revcomp(seq[start:end]) + seq[end:]
            truth.append(
                {"pos": orig, "end": orig + size, "svtype": "INV", "svlen": size}
            )
        elif svtype == "DUP":
            seq = seq[:end] + seq[start:end] + seq[end:]
            truth.append(
                {"pos": orig, "end": orig + size, "svtype": "DUP", "svlen": size}
            )
    return seq, truth


def write_vcf(
    path: Path, contigs: list[tuple[str, int]], records: dict[str, list[dict]]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##source=sv_spikein\n")
        fh.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">\n')
        fh.write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">\n')
        fh.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n')
        for name, length in contigs:
            fh.write(f"##contig=<ID={name},length={length}>\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        n = 0
        for name, _ in contigs:
            for rec in sorted(records.get(name, []), key=lambda r: r["pos"]):
                n += 1
                info = f"SVTYPE={rec['svtype']};SVLEN={rec['svlen']};END={rec['end']}"
                alt = f"<{rec['svtype']}>"
                fh.write(
                    f"{name}\t{rec['pos']}\tsv{n}\tN\t{alt}\t60\tPASS\t{info}\n"
                )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True, help="original reference FASTA")
    ap.add_argument("--variant-ref", required=True, help="donor reference FASTA")
    ap.add_argument("--truth-vcf", required=True, help="small-variant truth VCF")
    ap.add_argument("--output-ref", required=True, help="SV donor FASTA to write")
    ap.add_argument("--output-vcf", required=True, help="SV truth VCF to write")
    ap.add_argument("--seed", type=int, default=20260922)
    ap.add_argument(
        "--events-per-combo",
        type=int,
        default=EVENTS_PER_COMBO,
        help="events per (type, size) combination; 3 -> 39 events, 8 -> 104",
    )
    ap.add_argument(
        "--tag",
        default="",
        help="suffix inserted into output names, e.g. _s11",
    )
    ap.add_argument(
        "--force", action="store_true", help="overwrite outputs if they exist"
    )
    args = ap.parse_args()

    out_ref = Path(args.output_ref)
    out_vcf = Path(args.output_vcf)
    if args.tag:
        out_ref = out_ref.with_name(
            out_ref.stem + args.tag + out_ref.suffix
        )
        out_vcf = out_vcf.with_name(out_vcf.stem + args.tag + out_vcf.suffix)
    if (out_ref.exists() or out_vcf.exists()) and not args.force:
        print(f"outputs already exist: {out_ref} {out_vcf} (use --force to rewrite)")
        return

    rng = random.Random(args.seed)
    variants = parse_variants(Path(args.truth_vcf))
    offset_map = build_offset_map(variants)

    ref_seqs = read_fasta(Path(args.ref))
    variant_seqs = read_fasta(Path(args.variant_ref))
    lengths = {name: len(seq) for name, seq in ref_seqs}
    events = plan_events(lengths, rng, args.events_per_combo)

    donor: list[tuple[str, str]] = []
    truth: dict[str, list[dict]] = {}
    for name, seq in variant_seqs:
        entries = offset_map.get(name, [])
        contig_events = []
        for ev in events.get(name, []):
            vpos = to_variant_coord(entries, ev["pos"])
            if vpos < 1 or vpos + ev["size"] >= len(seq):
                continue
            contig_events.append({**ev, "pos": vpos})
        new_seq, recs = inject(seq, contig_events, rng)
        donor.append((name, new_seq))
        truth.setdefault(name, []).extend(recs)

    write_fasta(out_ref, donor)
    write_vcf(out_vcf, [(name, len(seq)) for name, seq in ref_seqs], truth)

    per_type: dict[str, int] = {}
    for recs in truth.values():
        for rec in recs:
            per_type[rec["svtype"]] = per_type.get(rec["svtype"], 0) + 1
    print(f"donor reference: {out_ref}")
    print(f"truth VCF: {out_vcf}")
    print(f"events: {sum(per_type.values())} {per_type}")


if __name__ == "__main__":
    main()
