#!/usr/bin/env python3
"""Dependency-free static checks for the reproducible release surface."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def run(command: list[str], cwd: Path) -> tuple[bool, str]:
    result = subprocess.run(
        command,
        cwd=cwd,
        capture_output=True,
        text=True,
    )
    output = (result.stdout + result.stderr).strip()
    return result.returncode == 0, output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--repo-root",
        default=str(Path(__file__).resolve().parents[1]),
    )
    args = parser.parse_args()
    root = Path(args.repo_root).resolve()

    failures: list[str] = []
    python_files = sorted(
        path
        for directory in ("src", "scripts", "tests")
        for path in (root / directory).rglob("*.py")
    )
    shell_files = sorted(
        path
        for directory in ("scripts",)
        for path in (root / directory).rglob("*.sh")
    )
    shell_files.append(root / "reproduce.sh")

    for path in python_files:
        ok, output = run(
            [sys.executable, "-m", "py_compile", str(path)],
            root,
        )
        if not ok:
            failures.append(f"python syntax: {path.relative_to(root)}\n{output}")

    for path in shell_files:
        if not path.exists():
            failures.append(f"missing shell entry point: {path}")
            continue
        ok, output = run(["bash", "-n", str(path)], root)
        if not ok:
            failures.append(f"shell syntax: {path.relative_to(root)}\n{output}")

    required = (
        "README.md",
        "REPRODUCIBILITY.md",
        "reproduce.sh",
        "requirements.txt",
        "environment.yml",
        "environment-tools.yml",
        "config/species_config_v2.json",
        "config/data_sources.json",
    )
    for relative in required:
        if not (root / relative).exists():
            failures.append(f"missing required release file: {relative}")

    if failures:
        print("\n\n".join(failures))
        raise SystemExit(1)
    print(
        f"STATIC_AUDIT_OK python_files={len(python_files)} "
        f"shell_files={len(shell_files)}"
    )


if __name__ == "__main__":
    main()
