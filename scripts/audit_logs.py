#!/usr/bin/env python3
"""Summarise failures and integrity signals in recent pipeline logs."""

from __future__ import annotations

import argparse
import re
import time
from pathlib import Path


PATTERNS = {
    "traceback": re.compile(r"Traceback \(most recent call last\)"),
    "crc32": re.compile(r"CRC32 checksum mismatch"),
    "missing_file": re.compile(r"FileNotFoundError|No such file or directory"),
    "failure": re.compile(r"\bERROR\b|\bFAILED\b|failed=", re.IGNORECASE),
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--days", type=float, default=3.0)
    parser.add_argument("--max-bytes", type=int, default=20 * 1024 * 1024)
    args = parser.parse_args()

    root = Path(args.project_dir).expanduser().resolve()
    cutoff = time.time() - args.days * 86400
    files = [
        path
        for path in (root / "logs").rglob("*")
        if path.is_file()
        and path.stat().st_mtime >= cutoff
        and path.stat().st_size <= args.max_bytes
    ]

    findings = []
    for path in sorted(files, key=lambda item: item.stat().st_mtime, reverse=True):
        text = path.read_text(encoding="utf-8", errors="replace")
        for category, pattern in PATTERNS.items():
            matches = [
                (line_number, line.strip())
                for line_number, line in enumerate(text.splitlines(), 1)
                if pattern.search(line)
            ]
            if matches:
                line_number, line = matches[-1]
                findings.append(
                    (path.relative_to(root), category, line_number, line[:240])
                )

    lines = [
        "# Recent log audit",
        "",
        f"project={root}",
        f"window_days={args.days} files_scanned={len(files)}",
        f"finding_groups={len(findings)}",
        "",
    ]
    for path, category, line_number, line in findings:
        lines.append(f"- `{path}:{line_number}` [{category}] {line}")
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"LOG_AUDIT={output} files={len(files)} findings={len(findings)}")
    for item in findings[:100]:
        print(f"{item[0]}:{item[2]} [{item[1]}] {item[3]}")


if __name__ == "__main__":
    main()
