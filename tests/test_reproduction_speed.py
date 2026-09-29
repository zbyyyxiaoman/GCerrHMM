"""Tests for deterministic quick-mode settings and release packaging inputs."""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from common_io import sha256_file  # noqa: E402
from package_github_release import CORE_SCRIPTS  # noqa: E402
from reproduce_gc_improvement import (  # noqa: E402
    allocate_threads,
    resolve_run_settings,
    resolve_threads,
)
from prepare_demo_data import (  # noqa: E402
    archive_corrupt,
    gzip_complete,
    parse_size,
    select_sort_threads,
)


def test_quick_settings_are_small_and_deterministic():
    settings = resolve_run_settings(quick=True)
    assert settings["coverage"] == 1.0
    assert settings["max_reads"] == 1000
    assert settings["evaluation_reads"] == 500
    assert settings["kmer_reads"] == 50
    assert settings["kmer_sample_size"] == 2000


def test_explicit_values_override_quick_defaults():
    settings = resolve_run_settings(
        quick=True,
        coverage=2.5,
        max_reads=17,
        evaluation_reads=11,
        kmer_reads=7,
        kmer_sample_size=101,
    )
    assert settings == {
        "coverage": 2.5,
        "max_reads": 17,
        "evaluation_reads": 11,
        "kmer_reads": 7,
        "kmer_sample_size": 101,
    }


def test_thread_budget_is_capped_and_allocated():
    assert resolve_threads(64) == 16
    assert resolve_threads(0) == 1
    assert allocate_threads(16, 3) == [6, 5, 5]
    assert allocate_threads(2, 3) == [1, 1]


def test_sort_memory_and_threads_are_bounded_by_available_memory():
    assert parse_size("1G") == 1 << 30
    assert parse_size("512M") == 512 << 20
    assert select_sort_threads(
        16,
        1 << 30,
        available_bytes=8 << 30,
    ) == 4
    assert select_sort_threads(
        16,
        1 << 30,
        explicit_sort_threads=2,
        available_bytes=8 << 30,
    ) == 2


def test_corrupt_gzip_is_rejected_and_archived():
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "broken.fastq.gz"
        path.write_bytes(b"\x1f\x8b\x08\x00truncated")
        assert not gzip_complete(path)
        archived = archive_corrupt(path)
        assert not path.exists()
        assert archived.exists()


def test_release_script_manifest_matches_the_checkout():
    missing = [
        filename
        for filename in CORE_SCRIPTS
        if not (ROOT / "scripts" / filename).exists()
    ]
    assert not missing


def test_tiny_fixture_is_present():
    fixture = ROOT / "tests" / "fixtures" / "tiny_demo" / "data"
    for relative in (
        "references/Ecoli_ref.fa",
        "real_reads_verified/Ecoli_ont.fastq.gz",
        "real_reads_verified/Ecoli_ont_aligned.bam",
        "real_reads_verified/Ecoli_ont_aligned.bam.bai",
    ):
        path = fixture / relative
        assert path.exists() and path.stat().st_size > 0


def test_sha256_file_reads_in_blocks():
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "payload.bin"
        path.write_bytes(b"abc")
        assert sha256_file(path, block_size=1) == (
            "ba7816bf8f01cfea414140de5dae2223"
            "b00361a396177a9cb410ff61f20015ad"
        )
