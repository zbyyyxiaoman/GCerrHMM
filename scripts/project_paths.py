"""Path helpers for repositories with separate code and result roots."""

from __future__ import annotations

from pathlib import Path


def project_file(project_dir: Path, relative: str) -> Path:
    candidates = (
        project_dir / relative,
        project_dir / "code" / relative,
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]
