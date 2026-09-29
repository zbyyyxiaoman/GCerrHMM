"""Matching logic checks for the SV spike-in evaluator."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from sv_eval import match, prf  # noqa: E402


def rec(chrom, pos, end, svtype, svlen):
    return {
        "chrom": chrom,
        "pos": pos,
        "end": end,
        "svtype": svtype,
        "svlen": svlen,
    }


def test_deletion_overlap_is_true_positive():
    truth = [rec("chr1", 1000, 6000, "DEL", 5000)]
    calls = [rec("chr1", 1010, 5990, "DEL", 4980)]
    tp, fp, fn = match(truth, calls)
    assert (len(tp), len(fp), len(fn)) == (1, 0, 0)


def test_type_mismatch_is_not_a_match():
    truth = [rec("chr1", 1000, 6000, "DEL", 5000)]
    calls = [rec("chr1", 1000, 6000, "DUP", 5000)]
    tp, fp, fn = match(truth, calls)
    assert (len(tp), len(fp), len(fn)) == (0, 1, 1)


def test_small_overlap_is_not_a_match():
    truth = [rec("chr1", 1000, 6000, "DEL", 5000)]
    calls = [rec("chr1", 5800, 9000, "DEL", 3200)]
    tp, fp, fn = match(truth, calls)
    assert (len(tp), len(fp), len(fn)) == (0, 1, 1)


def test_insertion_needs_position_and_length_agreement():
    truth = [rec("chr1", 5000, 5000, "INS", 1000)]
    near = [rec("chr1", 5100, 5100, "INS", 950)]
    far = [rec("chr1", 9000, 9000, "INS", 1000)]
    short = [rec("chr1", 5000, 5000, "INS", 100)]
    assert match(truth, near)[0]
    assert not match(truth, far)[0]
    assert not match(truth, short)[0]


def test_each_call_used_once():
    truth = [
        rec("chr1", 1000, 6000, "DEL", 5000),
        rec("chr1", 1100, 6100, "DEL", 5000),
    ]
    calls = [rec("chr1", 1050, 6050, "DEL", 5000)]
    tp, fp, fn = match(truth, calls)
    assert (len(tp), len(fp), len(fn)) == (1, 0, 1)


def test_prf_handles_empty_sets():
    assert prf(0, 0, 0)["f1"] == 0.0
    assert prf(5, 0, 0)["f1"] == 1.0
