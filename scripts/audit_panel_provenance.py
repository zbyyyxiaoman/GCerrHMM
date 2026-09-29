#!/usr/bin/env python3
"""Fail loudly if the delta-to-real panel and exported tables diverge.

The project has two separate chr21 data products:

* coverage-matched 30x files: ``*_Hsapiens_chr21_30x.json``
* the broader framework panel: ``*_alignment.json`` and ``*_fair_*.json``

The pre-registered decision rule names the first. This audit checks that the
decision, its standalone figure table, and the exported framework CSV rows all
resolve to the same 30x values. It is intentionally independent of the
producer scripts so a silent source-selection change cannot pass unnoticed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


SPECIES = "Hsapiens_chr21"
SUBSPECIES = "30x"
TOOLS = [
    ("errhmm", "GCerrHMM"),
    ("nanosim", "NanoSim"),
    ("pbsim", "PBSim3-sample"),
    ("badread", "badread"),
]
REAL = "real_ont30"


def read_json(path: Path) -> dict:
    if not path.exists() or path.stat().st_size == 0:
        raise SystemExit(f"missing provenance input: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path: Path) -> list[dict]:
    if not path.exists() or path.stat().st_size == 0:
        raise SystemExit(f"missing derived table: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def close(left, right, tolerance: float = 1e-9) -> bool:
    try:
        return math.isclose(float(left), float(right), rel_tol=0, abs_tol=tolerance)
    except (TypeError, ValueError):
        return False


def average_rank(values: dict[str, float]) -> dict[str, float]:
    ordered = sorted(values.items(), key=lambda item: item[1])
    result: dict[str, float] = {}
    index = 0
    while index < len(ordered):
        end = index + 1
        while end < len(ordered) and ordered[end][1] == ordered[index][1]:
            end += 1
        rank = (index + 1 + end) / 2
        for key, _ in ordered[index:end]:
            result[key] = rank
        index = end
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", required=True)
    args = parser.parse_args()

    project = Path(args.project_dir).expanduser()
    tables = project / "results" / "framework" / "tables"
    stats = project / "results" / "framework" / "stats"

    panel: dict[str, dict[str, dict]] = {}
    panel_paths: dict[str, dict[str, Path]] = {}
    for key, label in TOOLS:
        paths = {
            "mapping": tables / f"mapping_{key}_{SPECIES}_{SUBSPECIES}.json",
            "variant": tables / f"variant_{key}_{SPECIES}_{SUBSPECIES}.json",
            "assembly": tables / f"assembly_flye_{key}_{SPECIES}_{SUBSPECIES}.json",
        }
        panel[key] = {kind: read_json(path) for kind, path in paths.items()}
        panel_paths[key] = paths

    real = {
        "mapping": read_json(tables / f"mapping_{REAL}_{SPECIES}.json"),
        "variant": read_json(tables / f"variant_{REAL}_{SPECIES}.json"),
        "assembly": read_json(tables / f"assembly_flye_{REAL}_{SPECIES}.json"),
    }

    decision_rows = {
        row["tool_key"]: row
        for row in read_csv(stats / "delta_to_real_decision.csv")
    }
    if set(decision_rows) != {key for key, _ in TOOLS}:
        raise SystemExit("delta_to_real_decision.csv has an unexpected tool set")

    for key, label in TOOLS:
        row = decision_rows[key]
        expected = {
            "alignment_identity": panel[key]["mapping"]["base_identity"],
            "variant_f1": panel[key]["variant"]["f1_score"],
            "assembly_identity": panel[key]["assembly"]["reference_identity"],
            "assembly_n50": panel[key]["assembly"]["n50"],
            "assembly_contigs": panel[key]["assembly"]["num_contigs"],
        }
        for field, expected_value in expected.items():
            if not close(row.get(field), expected_value):
                raise SystemExit(
                    f"decision mismatch for {label} {field}: "
                    f"{row.get(field)!r} != {expected_value!r}"
                )

    r3_distances = {
        key: abs(
            panel[key]["mapping"]["base_identity"]
            - real["mapping"]["base_identity"]
        )
        for key, _ in TOOLS
    }
    r3_ranks = average_rank(r3_distances)
    identity_ranks = average_rank({
        key: abs(
            panel[key]["assembly"]["reference_identity"]
            - real["assembly"]["reference_identity"]
        )
        for key, _ in TOOLS
    })
    n50_ranks = average_rank({
        key: abs(math.log2(
            panel[key]["assembly"]["n50"] / real["assembly"]["n50"]
        ))
        for key, _ in TOOLS
    })
    borda = {
        key: (identity_ranks[key] + n50_ranks[key]) / 2
        for key, _ in TOOLS
    }
    r5_ranks = average_rank(borda)

    for key, label in TOOLS:
        if not close(decision_rows[key]["rank_R3"], r3_ranks[key], 1e-6):
            raise SystemExit(f"R3 rank mismatch for {label}")
        if not close(decision_rows[key]["rank_R5"], r5_ranks[key], 1e-6):
            raise SystemExit(f"R5 rank mismatch for {label}")

    errhmm = decision_rows["errhmm"]
    if int(errhmm["top2_layers"]) != 2:
        raise SystemExit(
            f"expected GCerrHMM top-2 in two anchored layers, got "
            f"{errhmm['top2_layers']}"
        )
    if not close(r3_ranks["errhmm"], 2.0, 1e-6):
        raise SystemExit(f"expected R3 rank 2, got {r3_ranks['errhmm']}")
    if not close(r5_ranks["errhmm"], 1.0, 1e-6):
        raise SystemExit(f"expected R5 rank 1, got {r5_ranks['errhmm']}")

    panel_rows = {
        row["tool_key"]: row
        for row in read_csv(stats / "delta_to_real_panel_30x.csv")
        if row["tool_key"] != REAL
    }
    for key, label in TOOLS:
        if not close(panel_rows[key]["identity"], panel[key]["mapping"]["base_identity"]):
            raise SystemExit(f"standalone figure identity mismatch for {label}")
        if not close(panel_rows[key]["asm_n50"], panel[key]["assembly"]["n50"]):
            raise SystemExit(f"standalone figure N50 mismatch for {label}")

    exported = {
        "mapping": {
            (row["species"], row["tool_key"]): row
            for row in read_csv(stats / "mapping.csv")
        },
        "variant": {
            (row["species"], row["tool_key"]): row
            for row in read_csv(stats / "variant.csv")
        },
        "assembly": {
            (row["species"], row["tool_key"]): row
            for row in read_csv(stats / "assembly.csv")
        },
    }
    for key, label in TOOLS:
        checks = (
            ("mapping", "base_identity_mean", panel[key]["mapping"]["base_identity"]),
            ("variant", "f1_score_mean", panel[key]["variant"]["f1_score"]),
            (
                "assembly",
                "reference_identity_mean",
                panel[key]["assembly"]["reference_identity"],
            ),
            ("assembly", "n50_mean", panel[key]["assembly"]["n50"]),
            ("assembly", "num_contigs_mean", panel[key]["assembly"]["num_contigs"]),
        )
        for kind, field, expected_value in checks:
            row = exported[kind].get((SPECIES, key))
            if row is None or not close(row.get(field), expected_value):
                raise SystemExit(
                    f"exported {kind}.csv mismatch for {label} {field}: "
                    f"{None if row is None else row.get(field)!r} != "
                    f"{expected_value!r}"
                )

    print("PANEL_PROVENANCE_OK")
    print(f"tools={len(TOOLS)} species={SPECIES} coverage={SUBSPECIES}")
    print(f"GCerrHMM R3 rank={r3_ranks['errhmm']} R5 rank={r5_ranks['errhmm']}")
    for key, label in TOOLS:
        path = panel_paths[key]["mapping"]
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        print(f"{label} mapping_sha256={digest}")


if __name__ == "__main__":
    main()
