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
import os
import re
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

from common_io import md5_file


def parse_size(value: str) -> int:
    match = re.fullmatch(r"\s*([0-9]+(?:\.[0-9]+)?)\s*([KMGT]?)(?:B)?\s*", value, re.I)
    if not match:
        raise ValueError(f"invalid size: {value!r}")
    number = float(match.group(1))
    suffix = match.group(2).upper()
    scale = {"": 1, "K": 1 << 10, "M": 1 << 20, "G": 1 << 30, "T": 1 << 40}
    return int(number * scale[suffix])


def available_memory_bytes() -> int | None:
    meminfo = Path("/proc/meminfo")
    if meminfo.exists():
        for line in meminfo.read_text(encoding="ascii").splitlines():
            if line.startswith("MemAvailable:"):
                kib = int(line.split()[1])
                return kib * 1024
    try:
        pages = os.sysconf("SC_AVPHYS_PAGES")
        page_size = os.sysconf("SC_PAGE_SIZE")
        return int(pages) * int(page_size)
    except (AttributeError, OSError, ValueError):
        return None


def select_sort_threads(
    requested_threads: int,
    sort_memory_bytes: int,
    explicit_sort_threads: int | None = None,
    available_bytes: int | None = None,
) -> int:
    maximum = max(1, int(requested_threads) - 1)
    if explicit_sort_threads is not None:
        return max(1, min(int(explicit_sort_threads), maximum))
    available = (
        available_bytes
        if available_bytes is not None
        else available_memory_bytes()
    )
    if available is None:
        return min(maximum, 3)
    memory_budget = max(sort_memory_bytes, int(available * 0.5))
    return max(1, min(maximum, memory_budget // sort_memory_bytes))


def check_free_space(project: Path, required_bytes: int) -> None:
    free = shutil.disk_usage(project).free
    if free < required_bytes:
        raise SystemExit(
            "insufficient free disk space: "
            f"need {required_bytes / (1 << 30):.2f} GiB, "
            f"available {free / (1 << 30):.2f} GiB under {project}"
        )
    print(
        f"[space] free={free / (1 << 30):.2f} GiB "
        f"required={required_bytes / (1 << 30):.2f} GiB"
    )


def gzip_complete(path: Path) -> bool:
    try:
        with gzip.open(path, "rb") as handle:
            for _ in iter(lambda: handle.read(8 << 20), b""):
                pass
        return path.stat().st_size > 0
    except (OSError, EOFError):
        return False


def archive_corrupt(path: Path) -> Path:
    archive = path.with_name(path.name + ".corrupt")
    index = 1
    while archive.exists():
        archive = path.with_name(f"{path.name}.corrupt.{index}")
        index += 1
    path.replace(archive)
    print(f"[archive] corrupt input moved to {archive}")
    return archive


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
        expected_size = item.get("fastq_size")
        expected_md5 = item.get("fastq_md5")
        valid = (
            (expected_size is None or fastq.stat().st_size == expected_size)
            and (
                expected_md5 is None
                or md5_file(fastq).lower() == expected_md5.lower()
            )
            and gzip_complete(fastq)
        )
        if valid:
            print(f"[skip] verified FASTQ: {fastq}")
            return fastq
        archive_corrupt(fastq)
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
    sort_memory: str,
    sort_threads: int | None,
) -> Path:
    bam = (
        project
        / "data"
        / "real_reads_verified"
        / f"{species}_ont_aligned.bam"
    )
    for tool in ("minimap2", "samtools"):
        if shutil.which(tool) is None:
            raise SystemExit(
                f"missing {tool}; install the conda environment with "
                "`conda env create -f environment.yml`"
            )
    if bam.exists() and bam.stat().st_size > 0:
        check = subprocess.run(
            ["samtools", "quickcheck", "-v", str(bam)],
            capture_output=True,
            text=True,
        )
        if check.returncode == 0:
            print(f"[skip] verified BAM: {bam}")
            return bam
        archive_corrupt(bam)
        index = bam.with_suffix(bam.suffix + ".bai")
        if index.exists():
            archive_corrupt(index)

    threads = max(1, int(threads))
    sort_memory_bytes = parse_size(sort_memory)
    sort_threads = select_sort_threads(
        threads,
        sort_memory_bytes,
        explicit_sort_threads=sort_threads,
    )
    bam.parent.mkdir(parents=True, exist_ok=True)
    print(
        f"[align] minimap2 threads={threads}; "
        f"samtools sort threads={sort_threads} memory={sort_memory}"
    )
    mapper = subprocess.Popen(
        [
            "minimap2",
            "-ax",
            "map-ont",
            "-t",
            str(threads),
            str(reference),
            str(fastq),
        ],
        stdout=subprocess.PIPE,
    )
    sorter = subprocess.Popen(
        [
            "samtools",
            "sort",
            "-@",
            str(sort_threads),
            "-m",
            sort_memory,
            "-o",
            str(bam),
            "-",
        ],
        stdin=mapper.stdout,
    )
    mapper.stdout.close()
    sort_status = sorter.wait()
    map_status = mapper.wait()
    if map_status != 0 or sort_status != 0:
        raise SystemExit(
            f"alignment pipeline failed: minimap2={map_status}, "
            f"samtools_sort={sort_status}"
        )
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
    parser.add_argument("--sort-memory", default="1G")
    parser.add_argument("--sort-threads", type=int, default=None)
    parser.add_argument("--min-free-gb", type=float, default=None)
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
    if args.species in {"Mmusculus", "Hsapiens"}:
        raise SystemExit(
            f"{args.species} is not the paper panel. The paper uses "
            "M. musculus chr19 and HG002 H. sapiens chr21 30x; see "
            "docs/panel_support.md for the panel-specific preparation."
        )

    source_bytes = (
        source_config["species"][args.species].get("sra_size")
        or source_config["species"][args.species].get("fastq_size")
        or 1 << 30
    )
    required_bytes = (
        int(args.min_free_gb * (1 << 30))
        if args.min_free_gb is not None
        else int(source_bytes * 2 + (1 << 30))
    )
    check_free_space(project, required_bytes)

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
        args.sort_memory,
        args.sort_threads,
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
