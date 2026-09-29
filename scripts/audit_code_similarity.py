#!/usr/bin/env python3
"""Report copied functions, repeated helpers and near-duplicate Python files."""

from __future__ import annotations

import argparse
import ast
import difflib
import hashlib
import re
from collections import Counter, defaultdict
from pathlib import Path


def python_files(root: Path) -> list[Path]:
    return sorted(
        path
        for directory in ("src", "scripts", "tests")
        for path in (root / directory).rglob("*.py")
    )


def normalized_lines(path: Path) -> list[str]:
    lines: list[str] = []
    in_docstring = False
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.count('"""') % 2 == 1 or line.count("'''") % 2 == 1:
            in_docstring = not in_docstring
            continue
        if in_docstring:
            continue
        line = re.sub(r"\s+", " ", line)
        lines.append(line)
    return lines


def function_hash(node: ast.AST) -> str:
    payload = ast.dump(node, annotate_fields=False, include_attributes=False)
    return hashlib.sha256(payload.encode()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    parser.add_argument("--output")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    files = python_files(root)

    function_locations: dict[str, list[str]] = defaultdict(list)
    function_names: Counter[str] = Counter()
    file_lines: dict[Path, list[str]] = {}
    for path in files:
        relative = path.relative_to(root).as_posix()
        file_lines[path] = normalized_lines(path)
        try:
            tree = ast.parse(path.read_text(encoding="utf-8", errors="replace"))
        except SyntaxError:
            continue
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                function_names[node.name] += 1
                function_locations[function_hash(node)].append(
                    f"{relative}:{node.lineno}:{node.name}"
                )

    exact_duplicates = {
        digest: locations
        for digest, locations in function_locations.items()
        if len(locations) > 1
    }
    repeated_names = {
        name: count for name, count in function_names.items() if count > 1
    }

    near_duplicates: list[tuple[float, str, str]] = []
    paths = list(file_lines)
    for index, left in enumerate(paths):
        for right in paths[index + 1:]:
            ratio = difflib.SequenceMatcher(
                None,
                file_lines[left],
                file_lines[right],
                autojunk=False,
            ).ratio()
            if ratio >= 0.70:
                near_duplicates.append(
                    (
                        ratio,
                        left.relative_to(root).as_posix(),
                        right.relative_to(root).as_posix(),
                    )
                )
    near_duplicates.sort(reverse=True)

    line_windows: dict[str, list[str]] = defaultdict(list)
    for path, lines in file_lines.items():
        relative = path.relative_to(root).as_posix()
        for index in range(0, max(0, len(lines) - 5)):
            window = "\n".join(lines[index:index + 6])
            digest = hashlib.sha256(window.encode()).hexdigest()
            line_windows[digest].append(f"{relative}:{index + 1}")
    repeated_windows = sorted(
        (
            (locations, window)
            for window, locations in line_windows.items()
            if len(locations) > 1
        ),
        key=lambda item: len(item[0]),
        reverse=True,
    )

    lines = [
        "# Code similarity audit",
        "",
        f"Scanned Python files: {len(files)}",
        "",
        "## Exact duplicate function bodies",
        "",
    ]
    if exact_duplicates:
        for locations in exact_duplicates.values():
            lines.append("- " + "; ".join(locations))
    else:
        lines.append("- None")
    lines.extend(["", "## Repeated helper names", ""])
    if repeated_names:
        for name, count in sorted(repeated_names.items()):
            lines.append(f"- `{name}`: {count} definitions")
    else:
        lines.append("- None")
    lines.extend(["", "## Near-duplicate files (>= 0.70)", ""])
    if near_duplicates:
        for ratio, left, right in near_duplicates:
            lines.append(f"- {ratio:.3f}: `{left}` vs `{right}`")
    else:
        lines.append("- None")
    lines.extend(["", "## Repeated six-line code windows", ""])
    if repeated_windows:
        for locations, window in repeated_windows[:20]:
            preview = " / ".join(window.splitlines()[:2])[:180]
            lines.append(f"- {len(locations)} copies: {preview}")
    else:
        lines.append("- None")

    report = "\n".join(lines) + "\n"
    if args.output:
        Path(args.output).write_text(report, encoding="utf-8")
    print(report)


if __name__ == "__main__":
    main()
