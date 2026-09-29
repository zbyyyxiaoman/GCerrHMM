#!/usr/bin/env python3
"""Download all configured source files and reuse verified legacy copies."""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

from common_io import first_fastq_read_id, md5_file


def verified(path, size, expected_md5=None):
    if not path.exists() or path.stat().st_size != size:
        return False
    return expected_md5 is None or md5_file(path).lower() == expected_md5.lower()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-dir", default=".")
    parser.add_argument("--config", default="config/data_sources.json")
    parser.add_argument(
        "--species",
        action="append",
        default=[],
        help="Download only the named species; may be supplied more than once.",
    )
    parser.add_argument("--parts", type=int, default=8)
    parser.add_argument("--chunk-size", type=int, default=16 << 20)
    args = parser.parse_args()

    project = Path(args.project_dir).resolve()
    with open(project / args.config, encoding="utf-8") as handle:
        config = json.load(handle)

    species_items = config["species"]
    if args.species:
        unknown = sorted(set(args.species) - set(species_items))
        if unknown:
            raise SystemExit(
                "unknown species: " + ", ".join(unknown)
                + "; available: " + ", ".join(sorted(species_items))
            )
        species_items = {
            key: species_items[key] for key in args.species
        }

    downloader = project / "scripts" / "download_fastq_multipart.py"
    for species, item in species_items.items():
        fastq = project / config["download_root"] / f"{species}_ont.fastq.gz"
        if item.get("fastq_verification") == "sra_conversion":
            sra = project / config["sra_root"] / f"{item['accession']}.sra"
            if not verified(sra, item["sra_size"]):
                subprocess.run([
                    sys.executable, str(downloader),
                    "--url", item["sra_url"],
                    "--output", str(sra),
                    "--size", str(item["sra_size"]),
                    "--parts", str(args.parts),
                    "--chunk-size", str(args.chunk_size),
                ], check=True)
            continue

        if verified(fastq, item["fastq_size"], item["fastq_md5"]):
            print(f"[SKIP] {fastq}")
            continue

        legacy = item.get("legacy_fastq")
        if item.get("fastq_verification") == "legacy_fastq" and legacy:
            legacy_path = project / legacy
            if (
                legacy_path.exists()
                and first_fastq_read_id(legacy_path).lstrip("@").startswith(
                    item["accession"]
                )
            ):
                if fastq.exists() or fastq.is_symlink():
                    fastq.unlink()
                os.link(legacy_path, fastq)
                print(f"[LINK] {fastq} -> {legacy_path}")
                continue

        if legacy:
            legacy_path = project / legacy
            if verified(legacy_path, item["fastq_size"], item["fastq_md5"]):
                if fastq.exists() or fastq.is_symlink():
                    fastq.unlink()
                os.link(legacy_path, fastq)
                print(f"[LINK] {fastq} -> {legacy_path}")
                continue

        subprocess.run([
            sys.executable, str(downloader),
            "--url", item["fastq_url"],
            "--output", str(fastq),
            "--size", str(item["fastq_size"]),
            "--md5", item["fastq_md5"],
            "--parts", str(args.parts),
            "--chunk-size", str(args.chunk_size),
        ], check=True)

    return 0


if __name__ == "__main__":
    sys.exit(main())
