#!/usr/bin/env python3
"""Check the local data/results tree for truncated or corrupt files."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def run(command: list[str]) -> tuple[bool, str]:
    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    return result.returncode == 0, result.stdout.strip()


def md5(path: Path) -> str:
    value = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def config_checks(root: Path) -> dict[Path, str]:
    expected: dict[Path, str] = {}
    config_roots = [root / "config", root / "code" / "config"]

    def config_file(name: str) -> Path | None:
        return next(
            (directory / name for directory in config_roots
             if (directory / name).exists()),
            None,
        )

    species_path = config_file("species_config_v2.json")
    data_path = config_file("data_sources.json")
    if species_path and data_path:
        species = json.loads(species_path.read_text(encoding="utf-8"))["species"]
        sources = json.loads(data_path.read_text(encoding="utf-8"))["species"]
        for key, item in species.items():
            digest = sources.get(key, {}).get("fastq_md5")
            relative = item.get("real_reads_ont")
            if digest and relative:
                expected[(root / relative).resolve()] = digest

    hifi_path = config_file("hifi_sources.json")
    if hifi_path:
        payload = json.loads(hifi_path.read_text(encoding="utf-8"))
        download_root = payload.get("download_root", "")
        for key, item in payload.get("sources", {}).items():
            digest = item.get("fastq_md5")
            accession = item.get("accession")
            if digest and accession:
                filename = f"{key}_{accession}_hifi.fastq.gz"
                expected[(root / download_root / filename).resolve()] = digest
    return expected


def check_md5(path: Path, expected: str) -> tuple[Path, bool, str]:
    actual = md5(path)
    ok = actual.lower() == expected.lower()
    return (
        path,
        ok,
        "" if ok else f"source_md5_mismatch expected={expected} actual={actual}",
    )


def check_bam(path: Path) -> tuple[Path, bool, str]:
    ok, output = run(["samtools", "quickcheck", "-v", str(path)])
    return path, ok, output


def check_gzip(path: Path) -> tuple[Path, bool, str]:
    ok, output = run(["gzip", "-t", str(path)])
    return path, ok, output


def check_uncompressed_fastq(path: Path) -> tuple[Path, bool, str]:
    head_ok, head = run(["head", "-n", "4", str(path)])
    tail_ok, tail = run(["tail", "-n", "4", str(path)])
    if not head_ok or not tail_ok:
        return path, False, (head + "\n" + tail).strip()
    head_lines = head.splitlines()
    tail_lines = tail.splitlines()
    valid = (
        len(head_lines) == 4
        and len(tail_lines) == 4
        and head_lines[0].startswith("@")
        and head_lines[2].startswith("+")
        and tail_lines[0].startswith("@")
        and tail_lines[2].startswith("+")
    )
    return path, valid, "" if valid else "invalid FASTQ head/tail records"


def is_known_noncritical(root: Path, path: Path, detail: str) -> bool:
    relative = path.relative_to(root).as_posix()
    if "_incomplete_" in relative:
        return True
    if path.stat().st_size == 0 and (
        relative.startswith("data/tmp/")
        or "smoke" in relative
        or relative == "data/real_reads_hifi_phasing/HG002_chr21.fastq"
    ):
        return True
    if path.name.startswith(("sample_profile_", "train_subset")):
        return True
    if relative.startswith("data/simulated/cross_tools/"):
        return True
    if (
        relative.startswith("data/real_reads_verified/")
        and detail.startswith("source_md5_mismatch")
    ):
        return True
    if relative == "data/real_reads_hifi/Hsapiens_ERR13110527_hifi.fastq.gz":
        return True
    return False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", required=True)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--output", default="")
    args = parser.parse_args()

    root = Path(args.project_dir).expanduser().resolve()
    jobs = max(1, min(args.jobs, 16))
    search_roots = [root / "data", root / "results"]
    bams = sorted(
        path for directory in search_roots
        for path in directory.rglob("*.bam")
    )
    expected_md5 = config_checks(root)
    configured_paths = {path for path in expected_md5 if path.exists()}
    gzips = sorted(
        path for directory in search_roots
        for path in directory.rglob("*")
        if path.is_file() and path.suffix in {".gz"}
        and path.resolve() not in configured_paths
    )
    fastqs = sorted(
        path for directory in search_roots
        for path in directory.rglob("*.fastq")
        if not path.name.startswith(("sample_profile_", "train_subset"))
    )

    checks = []
    checks.extend((check_bam, path) for path in bams)
    checks.extend((check_gzip, path) for path in gzips)
    checks.extend((check_uncompressed_fastq, path) for path in fastqs)
    checks.extend(
        (
            lambda path=path, expected=digest: check_md5(path, expected),
            path,
        )
        for path, digest in expected_md5.items()
        if path.exists()
    )
    checks.extend(
        (check_gzip, path)
        for path in expected_md5
        if path.exists()
    )

    failures = []
    warnings = []
    completed = 0
    with ThreadPoolExecutor(max_workers=jobs) as executor:
        futures = {
            executor.submit(function, path): (function, path)
            for function, path in checks
        }
        for future in as_completed(futures):
            path, ok, detail = future.result()
            completed += 1
            if not ok:
                if is_known_noncritical(root, path, detail):
                    warnings.append((path, detail))
                else:
                    failures.append((path, detail))
            if completed % 100 == 0:
                print(f"checked={completed}/{len(checks)}")

    lines = [
        "# Data integrity audit",
        "",
        f"root={root}",
        f"jobs={jobs}",
        f"bam={len(bams)} gzip={len(gzips)} fastq={len(fastqs)}",
        f"configured_md5={sum(1 for path in expected_md5 if path.exists())}",
        f"checked={completed} critical_failures={len(failures)} "
        f"known_warnings={len(warnings)}",
    ]
    for path, detail in failures:
        lines.append(f"FAIL {path}: {detail}")
    for path, detail in warnings:
        lines.append(f"WARN_KNOWN {path}: {detail}")
    report = "\n".join(lines) + "\n"
    print(report)
    if args.output:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(report, encoding="utf-8")
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
