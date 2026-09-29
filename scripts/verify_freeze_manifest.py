#!/usr/bin/env python3
"""Verify checksum bullets recorded in docs/data_freeze_manifest.md."""

import argparse
import hashlib
import re
from pathlib import Path


CHECKSUM_PATTERNS = (
    re.compile(r"^- `([0-9a-fA-F]{64})  (.+)`$"),
    re.compile(r"^([0-9a-fA-F]{64})  (.+)$"),
)


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def manifest_entries(manifest: Path):
    """Yield the latest hash for each relative artifact path."""
    entries: dict[str, tuple[str, str]] = {}
    for line in manifest.read_text(encoding="utf-8").splitlines():
        match = next(
            (pattern.match(line.strip()) for pattern in CHECKSUM_PATTERNS
             if pattern.match(line.strip())),
            None,
        )
        if match is None:
            continue
        expected, raw_path = match.groups()
        entries[raw_path] = (expected.lower(), raw_path)
    return entries.values()


def resolve_artifact(root: Path, raw_path: str) -> tuple[Path, str]:
    normalized = raw_path.replace("\\", "/")
    candidate = Path(normalized)
    if candidate.is_absolute():
        for marker in (
            f"/{root.name}/",
            "/project/",
            "/errhmm_project/",
        ):
            if marker in normalized:
                relative = normalized.split(marker, 1)[1]
                return root / relative, relative
        return candidate, normalized.lstrip("/")
    return root / candidate, normalized


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", required=True)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()
    root = Path(args.project_dir)
    manifest = root / "docs" / "data_freeze_manifest.md"
    checked = 0
    failed = 0
    missing = 0
    for expected, raw_path in manifest_entries(manifest):
        path, relative = resolve_artifact(root, raw_path)
        if not path.exists():
            missing += 1
            print(f"MISSING {relative}")
            continue
        actual = digest(path)
        checked += 1
        if actual != expected.lower():
            failed += 1
            print(f"FAIL {relative}: expected={expected.lower()} actual={actual}")
        else:
            print(f"OK {relative}")
    print(
        f"FREEZE_MANIFEST checked={checked} failed={failed} missing={missing}"
    )
    if failed or (args.strict and missing):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
