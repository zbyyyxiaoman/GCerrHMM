#!/usr/bin/env python3
"""Print the block span (and first-record uoffset) of one reference in a BAI.

The BAI pseudo-bin (37450) stores, for every reference, the virtual offset of
the first and last record. Splitting that into BGZF block offset and in-block
offset lets a ranged fetch pull exactly the bytes that hold a chromosome's
records without downloading or repairing the whole file.

usage: bai_ref_span.py --bam a.bam --index a.bam.bai --ref chr21
       -> "start end uoffset"
"""

from __future__ import annotations

import argparse
import struct
import subprocess


def ref_names(bam_path: str) -> list[str]:
    header = subprocess.run(
        ["samtools", "view", "-H", bam_path],
        capture_output=True, text=True, check=True,
    ).stdout
    names = []
    for line in header.splitlines():
        if line.startswith("@SQ"):
            for field in line.split("\t"):
                if field.startswith("SN:"):
                    names.append(field[3:])
    return names


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam", required=True)
    parser.add_argument("--index", required=True)
    parser.add_argument("--ref", required=True)
    args = parser.parse_args()

    names = ref_names(args.bam)
    data = open(args.index, "rb").read()
    if data[:4] != b"BAI\x01":
        raise SystemExit(f"not a BAI file: {data[:4]!r}")

    off = 4
    (n_ref,) = struct.unpack_from("<i", data, off)
    off += 4

    spans: dict[str, tuple[int, int]] = {}
    for rid in range(n_ref):
        (n_bin,) = struct.unpack_from("<i", data, off)
        off += 4
        pseudo = None
        for _ in range(n_bin):
            (bin_id,) = struct.unpack_from("<I", data, off)
            off += 4
            (n_chunk,) = struct.unpack_from("<i", data, off)
            off += 4
            chunks = []
            for _ in range(n_chunk):
                beg, end = struct.unpack_from("<QQ", data, off)
                off += 16
                chunks.append((beg, end))
            if bin_id == 37450:
                pseudo = chunks
        (n_intv,) = struct.unpack_from("<i", data, off)
        off += 4 + 8 * n_intv

        name = names[rid] if rid < len(names) else f"ref{rid}"
        if pseudo:
            spans[name] = (pseudo[0][0], pseudo[-1][1])

    if args.ref not in spans:
        raise SystemExit(f"no pseudo-bin records for {args.ref}")

    # Coordinate-sorted BAM: the span must stop where the next reference with
    # records begins, otherwise the final BGZF blocks are cut mid-chromosome.
    ordered = [rid for rid in range(n_ref)
               if (names[rid] if rid < len(names) else f"ref{rid}") in spans]
    rid_of_ref = ordered[[names[r] if r < len(names) else f"ref{r}"
                          for r in ordered].index(args.ref)]
    later = [r for r in ordered if r > rid_of_ref]
    if later:
        nxt = names[later[0]] if later[0] < len(names) else f"ref{later[0]}"
        end_block = spans[nxt][0] >> 16
    else:
        end_block = spans[args.ref][1] >> 16

    first_vo = spans[args.ref][0]
    print(f"{first_vo >> 16} {end_block} {first_vo & 0xFFFF}")


if __name__ == "__main__":
    main()
