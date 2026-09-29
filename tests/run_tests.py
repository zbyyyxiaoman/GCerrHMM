#!/usr/bin/env python3
"""Run plain-function tests without requiring pytest."""

from __future__ import annotations

import importlib.util
import inspect
import sys
import traceback
from pathlib import Path


def load_module(path: Path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    failures = []
    passed = 0
    for path in sorted((root / "tests").glob("test_*.py")):
        module = load_module(path)
        for name, value in sorted(vars(module).items()):
            if not name.startswith("test_") or not inspect.isfunction(value):
                continue
            try:
                value()
                passed += 1
            except Exception:
                failures.append(
                    f"{path.relative_to(root)}::{name}\n{traceback.format_exc()}"
                )
    if failures:
        print("\n\n".join(failures))
        raise SystemExit(1)
    print(f"TESTS_OK passed={passed}")


if __name__ == "__main__":
    main()
