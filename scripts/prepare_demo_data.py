#!/usr/bin/env python3
"""Prepare the public inputs required by ``reproduce.sh --quick-demo``.

This script intentionally downloads only one species. It fetches the NCBI
reference, delegates FASTQ acquisition to ``download_verified_sources.py``,
and aligns the reads with minimap2/samtools.
"""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path


def download(url: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".part")
    print(f"[download] {url}")
    print(f"[download] -> {destination}")
    with urllib.request.urlopen(url, timeout=120) as source:
        with temporary.open("wb") as target:
            shutil.copyfileobj(source, target, length=8 << 20)
    temporary.replace(destination)


def prepare_reference(project: Path, config: dict, species: str) -> Path:
    item = config["species"][species]
    reference = project / item["ref_fasta"]
    if reference.exists() and reference.stat().st_size > 0:
        print(f"[skip] reference exists: {reference}")
        return reference
    if "ref_url" not in item:
        raise SystemExit(f"no reference URL for {species}")
    archive = reference.with_suffix(reference.suffix + ".gz")
    if not archive.exists():
        download(item["ref_url"], archive)
    reference.parent.mkdir(parents=True, exist_ok=True)
    print(f"[gunzip] {archive} -> {reference}")
    with gzip.open(archive, "rb") as source:
        with reference.open("wb") as target:
            shutil.copyfileobj(source, target, length=8 << 20)
    return reference


def prepare_fastq(
    project: Path,
    config: dict,
    species: str,
    parts: int,
    chunk_size: int,
    sra_tool: str | None,
    sra_threads: int,
) -> Path:
    item = config["species"][species]
    fastq = project / config["download_root"] / f"{species}_ont.fastq.gz"
    if fastq.exists() and fastq.stat().st_size > 0:
        print(f"[skip] FASTQ exists: {fastq}")
        return fastq
    command = [
        sys.executable,
        str(project / "scripts" / "download_verified_sources.py"),
        "--project-dir",
        str(project),
        "--species",
        species,
        "--parts",
        str(parts),
        "--chunk-size",
        str(chunk_size),
    ]
    if sra_tool:
        command += ["--sra-tool", sra_tool]
    command += ["--sra-threads", str(sra_threads)]
    subprocess.run(command, check=True)
    if not fastq.exists() or fastq.stat().st_size == 0:
        raise SystemExit(f"FASTQ was not created: {fastq}")
    return fastq


def prepare_alignment(
    project: Path,
    species: str,
    reference: Path,
    fastq: Path,
    threads: int,
) -> Path:
    bam = (
        project
        / "data"
        / "real_reads_verified"
        / f"{species}_ont_aligned.bam"
    )
    if bam.exists() and bam.stat().st_size > 0:
        print(f"[skip] BAM exists: {bam}")
        return bam
    for tool in ("minimap2", "samtools"):
        if shutil.which(tool) is None:
            raise SystemExit(
                f"missing {tool}; install the conda environment with "
                "`conda env create -f environment.yml`"
            )
    threads = max(1, int(threads))
    sort_threads = max(1, threads - 1)
    bam.parent.mkdir(parents=True, exist_ok=True)
    command = (
        f"minimap2 -ax map-ont -t {threads} "
        f"{reference} {fastq} | "
        f"samtools sort -@ {sort_threads} -m 1G -o {bam} -"
    )
    print("[align] " + command)
    subprocess.run(command, shell=True, check=True)
    subprocess.run(["samtools", "index", str(bam)], check=True)
    return bam


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-dir", default=".")
    parser.add_argument("--species", default="Ecoli")
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--parts", type=int, default=8)
    parser.add_argument("--chunk-size", type=int, default=16 << 20)
    parser.add_argument("--sra-tool", default=None)
    parser.add_argument("--sra-threads", type=int, default=8)
    args = parser.parse_args()

    project = Path(args.project_dir).resolve()
    with (project / "config" / "species_config_v2.json").open(
        encoding="utf-8"
    ) as handle:
        species_config = json.load(handle)
    with (project / "config" / "data_sources.json").open(
        encoding="utf-8"
    ) as handle:
        source_config = json.load(handle)
    if args.species not in species_config["species"]:
        raise SystemExit(
            f"unknown species {args.species!r}; available: "
            + ", ".join(sorted(species_config["species"]))
        )
    if args.species not in source_config["species"]:
        raise SystemExit(f"{args.species} has no configured ONT source")

    reference = prepare_reference(
        project,
        species_config,
        args.species,
    )
    fastq = prepare_fastq(
        project,
        source_config,
        args.species,
        args.parts,
        args.chunk_size,
        args.sra_tool,
        args.sra_threads,
    )
    bam = prepare_alignment(
        project,
        args.species,
        reference,
        fastq,
        args.threads,
    )
    print("DEMO_DATA_READY")
    print(f"REFERENCE={reference}")
    print(f"FASTQ={fastq}")
    print(f"BAM={bam}")
    print(
        "NEXT: PROJECT_DIR={} THREADS={} JOBS=3 "
        "bash reproduce.sh --quick-demo".format(project, args.threads)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
