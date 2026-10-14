#!/usr/bin/env python3
"""Install the frozen public result bundle without overwriting conflicts."""

from __future__ import annotations

import argparse
import hashlib
import shutil
import sys
import tempfile
import zipfile
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_manifest(bundle: Path) -> list[tuple[str, str]]:
    manifest = bundle / "MANIFEST.sha256"
    if not manifest.exists():
        raise SystemExit(f"bundle manifest missing: {manifest}")
    entries = []
    for line in manifest.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        digest, relative = line.split("  ", 1)
        entries.append((digest, relative))
    return entries


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-dir", required=True)
    parser.add_argument("--bundle-dir", required=True)
    parser.add_argument("--check-only", action="store_true")
    parser.add_argument("--strict", action="store_true")
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()

    project = Path(args.project_dir).resolve()
    bundle = Path(args.bundle_dir).resolve()
    temporary = None
    if not bundle.exists():
        archive = bundle.with_name(bundle.name + ".zip")
        if not archive.exists():
            raise SystemExit(
                f"result bundle missing: {bundle} or {archive}"
            )
        temporary = tempfile.TemporaryDirectory(prefix="gcerrhmm-results-")
        bundle = Path(temporary.name)
        with zipfile.ZipFile(archive) as handle:
            handle.extractall(bundle)
    entries = read_manifest(bundle)

    failures = []
    for expected, relative in entries:
        source = bundle / relative
        if not source.exists():
            failures.append(f"missing bundle file: {relative}")
        elif sha256_file(source) != expected:
            failures.append(f"bundle checksum mismatch: {relative}")
    if failures:
        raise SystemExit("\n".join(failures))
    if args.check_only:
        print(f"BUNDLED_RESULTS_OK files={len(entries)}")
        if temporary is not None:
            temporary.cleanup()
        return 0

    target_root = project / "results"
    copied = 0
    existing = 0
    conflicts = []
    for _, relative in entries:
        source = bundle / relative
        target = target_root / relative
        if target.exists():
            if sha256_file(target) == sha256_file(source):
                existing += 1
            else:
                conflicts.append(relative)
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, target)
        copied += 1

    if conflicts and args.strict:
        preview = "\n".join(f"- {item}" for item in conflicts[:20])
        raise SystemExit(
            "existing results conflict with the frozen bundle; refusing to "
            f"mix vintages ({len(conflicts)} conflicts):\n{preview}"
        )

    if conflicts and not args.quiet:
        preview = "\n".join(f"- {item}" for item in conflicts[:20])
        print(
            f"BUNDLED_RESULTS_WARNING existing_conflicts={len(conflicts)}; "
            "existing files were preserved",
            file=sys.stderr,
        )
        print(preview, file=sys.stderr)

    print(
        f"BUNDLED_RESULTS_INSTALLED copied={copied} "
        f"existing={existing} conflicts={len(conflicts)}"
    )
    if temporary is not None:
        temporary.cleanup()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
