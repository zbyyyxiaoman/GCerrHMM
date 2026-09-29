"""Tests for manifest parsing and provenance-rank helpers."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from audit_data_integrity import is_known_noncritical  # noqa: E402
from audit_panel_provenance import average_rank, close  # noqa: E402
from project_paths import project_file  # noqa: E402
from verify_freeze_manifest import manifest_entries, resolve_artifact  # noqa: E402


def test_average_rank_handles_ties():
    values = {"a": 1.0, "b": 1.0, "c": 3.0}
    assert average_rank(values) == {"a": 1.5, "b": 1.5, "c": 3.0}


def test_close_is_strict_about_missing_values():
    assert close(1.0, 1.0 + 1e-10)
    assert not close(None, 1.0)


def test_manifest_uses_latest_entry_for_each_path():
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "manifest.md"
        path.write_text(
            "- `"
            + "a" * 64
            + "  results/example.csv`\n"
            + "```\n"
            + "b" * 64
            + "  results/example.csv\n"
            + "```\n",
            encoding="utf-8",
        )
        entries = list(manifest_entries(path))
        assert len(entries) == 1
        assert entries[0][0] == "b" * 64


def test_resolve_artifact_maps_remote_project_paths():
    root = Path("/opt/repro/errhmm_project")
    path, relative = resolve_artifact(
        root,
        "/srv/repro/errhmm_project/results/example.json",
    )
    assert path == root / "results/example.json"
    assert relative == "results/example.json"


def test_source_md5_mismatch_is_a_known_warning():
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        path = root / "data/real_reads_verified/Ecoli_ont.fastq.gz"
        path.parent.mkdir(parents=True)
        path.write_bytes(b"converted")
        assert is_known_noncritical(
            root,
            path,
            "source_md5_mismatch expected=x actual=y",
        )


def test_incomplete_archive_is_noncritical():
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        path = root / "data/_incomplete_sample/file.fastq.gz"
        path.parent.mkdir(parents=True)
        path.write_bytes(b"incomplete")
        assert is_known_noncritical(root, path, "gzip failed")


def test_project_file_falls_back_to_code_subdirectory():
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        path = root / "code/docs/paper_figures/table1_species_panel.csv"
        path.parent.mkdir(parents=True)
        path.write_text("species\n", encoding="utf-8")
        assert project_file(
            root,
            "docs/paper_figures/table1_species_panel.csv",
        ) == path
