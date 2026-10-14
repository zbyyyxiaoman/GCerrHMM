#!/usr/bin/env python3
"""Build the portable reviewer-facing result bundle.

The original analysis tree is intentionally not redistributed wholesale. This
script selects the text result artifacts needed by the public reproduction
entry points, removes machine-local paths, and writes a SHA-256 manifest.
"""

from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path


INCLUDE_ROOTS = (
    "framework/tables",
    "framework/stats",
    "tables",
    "stats",
)

EXTRA_FILES = (
    "framework/hg002_chr21_20260922_204518/phasing/phasing_result.json",
)

TEXT_SUFFIXES = {
    ".csv",
    ".json",
    ".md",
    ".sha256",
    ".txt",
}

PORTABLE_REPLACEMENTS = (
    (re.compile(r"/mnt/c/Users/[^/]+/Desktop/[^/]+/?"), ""),
    (re.compile(r"[A-Za-z]:[\\/]Users[\\/][^\\/]+[\\/]Desktop[\\/][^\\/]+[\\/]?"), ""),
    (re.compile(r"/home/[^/\s]+/errhmm_project/?"), ""),
    (re.compile(r"/home/[^/\s]+/bmc_data/?"), "<external-data>/"),
    (re.compile(r"/home/[^/\s]+/miniconda3"), "$CONDA_PREFIX"),
    (re.compile(r"/home/[^/\s]+/apt-bcftools/root/usr/bin"), "<bcftools-prefix>"),
    (re.compile(r"\b[\w.-]+@10\.70\.5\.64\b"), "<remote-host>"),
    (re.compile(r"\b10\.70\.5\.64\b"), "<remote-host>"),
    (re.compile(r"\bbmc_project\b"), "project"),
)

PRIVATE_PATTERNS = (
    re.compile(r"/mnt/c/Users/27947(?:/|\\|$)"),
    re.compile(r"[A-Za-z]:[\\/]Users[\\/]27947(?:[\\/]|$)"),
    re.compile(r"/home/bobby(?:/|$)"),
    re.compile(r"/home/zby(?:/|$)"),
    re.compile(r"\b[\w.-]+@10\.70\.5\.64\b"),
    re.compile(r"\b10\.70\.5\.64\b"),
    re.compile(r"Desktop[\\/]project(?:[\\/]|$)"),
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_text(path: Path) -> str:
    payload = path.read_bytes()
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError:
        text = payload.decode("gb18030")
    return text.replace("\r\n", "\n").replace("\r", "\n")


def portable_text(text: str) -> str:
    for pattern, replacement in PORTABLE_REPLACEMENTS:
        text = pattern.sub(replacement, text)
    return text


def candidate_files(source: Path):
    for root_name in INCLUDE_ROOTS:
        root = source / root_name
        if not root.exists():
            continue
        for path in sorted(root.rglob("*")):
            if path.is_file() and path.suffix.lower() in TEXT_SUFFIXES:
                yield path.relative_to(source), path
    for relative_name in EXTRA_FILES:
        path = source / relative_name
        if path.exists() and path.is_file():
            yield Path(relative_name), path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-results", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    source = Path(args.source_results).resolve()
    output = Path(args.output).resolve()
    if output.exists():
        raise SystemExit(
            f"refusing to overwrite existing bundle: {output}"
        )
    output.mkdir(parents=True)

    files = list(candidate_files(source))
    if not files:
        raise SystemExit(f"no result files found under {source}")

    for relative, source_path in files:
        destination = output / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        text = portable_text(read_text(source_path))
        destination.write_text(text, encoding="utf-8")

    marker_findings = []
    manifest_lines = []
    for destination in sorted(output.rglob("*")):
        if not destination.is_file():
            continue
        relative = destination.relative_to(output)
        text = read_text(destination)
        for pattern in PRIVATE_PATTERNS:
            if pattern.search(text):
                marker_findings.append(
                    f"{relative.as_posix()}: {pattern.pattern}"
                )
                break
        manifest_lines.append(
            f"{sha256_file(destination)}  {relative.as_posix()}"
        )

    if marker_findings:
        raise SystemExit(
            "private path audit failed:\n- " + "\n- ".join(marker_findings)
        )

    (output / "MANIFEST.sha256").write_text(
        "\n".join(manifest_lines) + "\n",
        encoding="utf-8",
    )
    print(f"BUNDLE_FILES={len(manifest_lines)}")
    print(f"BUNDLE_DIR={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
