#!/usr/bin/env python3
"""Report non-portable absolute paths in code and scripts."""

from __future__ import annotations

import argparse
import re
from pathlib import Path


PATTERNS = {
    "wsl_local": re.compile(r"/mnt/c/Users/"),
    "windows_local": re.compile(r"C:\\Users\\", re.IGNORECASE),
    "windows_forward_slash": re.compile(r"C:/Users/", re.IGNORECASE),
    "linux_home": re.compile(r"/(?:home|Users)/[^/\s]+/"),
    "forbidden_shared_storage": re.compile(r"/data/[12]/"),
}

EXCLUDED = {
    "scripts/audit_hardcoded_paths.py",
    "scripts/package_github_release.py",
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    root = Path(args.project_dir).resolve()
    findings = []
    files = sorted(
        path
        for directory in ("src", "scripts", "tests")
        for path in (root / directory).rglob("*")
        if path.is_file() and path.suffix in {".py", ".sh", ".R", ".r"}
    )
    for path in files:
        relative = path.relative_to(root).as_posix()
        if relative in EXCLUDED:
            continue
        for line_number, line in enumerate(
            path.read_text(encoding="utf-8", errors="replace").splitlines(),
            1,
        ):
            for category, pattern in PATTERNS.items():
                if pattern.search(line):
                    findings.append(
                        (
                            relative,
                            line_number,
                            category,
                            line.strip()[:220],
                        )
                    )

    lines = [
        "# Hardcoded path audit",
        "",
        f"files_scanned={len(files)} findings={len(findings)}",
        "",
    ]
    for relative, line_number, category, line in findings:
        lines.append(f"- `{relative}:{line_number}` [{category}] {line}")
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"PATH_AUDIT={output} files={len(files)} findings={len(findings)}")
    if args.strict and findings:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
