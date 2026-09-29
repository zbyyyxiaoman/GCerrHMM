#!/usr/bin/env python3
"""Write a random subset of a FASTQ up to a target number of bases.

Used to bring simulated read sets to the same coverage as the real ONT
control (HG002 chr21, 6.62x), so an assembly comparison is not confounded by
the simulators having ~1.5x more coverage.

Reads are kept whole and chosen without replacement with a fixed seed, so the
subset is reproducible.

Usage:
  python3 subsample_fastq_to_bases.py --input in.fastq --target-bases 309228361 \
      --output out.fastq --seed 7
"""

from __future__ import annotations

import argparse
import gzip
import random
from pathlib import Path


def opener(path: Path, mode: str):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--target-bases", type=int, required=True)
    parser.add_argument("--seed", type=int, default=7)
    args = parser.parse_args()

    src = Path(args.input)
    dst = Path(args.output)
    if dst.exists() and dst.stat().st_size > 0:
        print(f"SKIP existing {dst}")
        return

    records = []
    total = 0
    with opener(src, "r") as handle:
        block = []
        for line in handle:
            block.append(line)
            if len(block) == 4:
                records.append(block)
                total += len(block[1].strip())
                block = []

    rng = random.Random(args.seed)
    rng.shuffle(records)
    kept = []
    running = 0
    for record in records:
        if running >= args.target_bases:
            break
        kept.append(record)
        running += len(record[1].strip())

    dst.parent.mkdir(parents=True, exist_ok=True)
    with opener(dst, "w") as handle:
        for record in kept:
            handle.writelines(record)

    print(
        f"{src.name}: {len(records)} reads / {total} bp -> "
        f"{len(kept)} reads / {running} bp ({100 * running / total:.1f}%)"
    )


if __name__ == "__main__":
    main()
