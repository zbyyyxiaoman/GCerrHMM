"""Small, dependency-free helpers shared by analysis and packaging scripts."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any


_MISSING = object()


def sha256_file(path: Path, block_size: int = 1024 * 1024) -> str:
    """Return the SHA-256 digest of *path* without loading it into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def md5_file(path: Path, block_size: int = 8 << 20) -> str:
    """Return the MD5 digest of *path* for public data-source verification."""
    digest = hashlib.md5()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(block_size), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text_maybe_gzip(path: Path, encoding: str = "utf-8"):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    return opener(path, "rt", encoding=encoding, errors="replace")


def first_fastq_read_id(path: Path) -> str:
    with open_text_maybe_gzip(path) as handle:
        return handle.readline().strip().split()[0]


def read_json(
    path: Path,
    default: Any = _MISSING,
    *,
    ignore_errors: bool = False,
) -> Any:
    """Read JSON, optionally returning *default* for missing/empty/bad input."""
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        if default is _MISSING:
            raise FileNotFoundError(path)
        return default
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        if ignore_errors and default is not _MISSING:
            return default
        raise


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    """Read a CSV file as a list of dictionaries."""
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_csv_table(path: Path) -> list[list[str]]:
    """Read a CSV file as a list of rows, including its header."""
    with Path(path).open(newline="", encoding="utf-8") as handle:
        return list(csv.reader(handle))
