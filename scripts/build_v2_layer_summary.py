#!/usr/bin/env python3
"""Export a compact v2 summary of every anchored evaluation layer."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


TOOLS = [
    ("errhmm", "GCerrHMM"),
    ("nanosim", "NanoSim"),
    ("pbsim", "PBSim3-sample"),
    ("badread", "badread"),
]


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", required=True)
    args = parser.parse_args()

    root = Path(args.project_dir).expanduser()
    stats = root / "results/framework/stats"
    tables = root / "results/framework/tables"

    panel: dict[str, dict[str, dict]] = {}
    for key, _ in TOOLS:
        panel[key] = {
            "mapping": read_json(
                tables / f"mapping_{key}_Hsapiens_chr21_30x.json"
            ),
            "variant": read_json(
                tables / f"variant_{key}_Hsapiens_chr21_30x.json"
            ),
            "assembly": read_json(
                tables / f"assembly_flye_{key}_Hsapiens_chr21_30x.json"
            ),
        }
    real_mapping = read_json(
        tables / "mapping_real_ont30_Hsapiens_chr21.json"
    )
    real_variant = read_json(
        tables / "variant_real_ont30_Hsapiens_chr21.json"
    )
    real_assembly = read_json(
        tables / "assembly_flye_real_ont30_Hsapiens_chr21.json"
    )

    rows = []
    for layer, metric, real_value, source_name, value_field, note in (
        (
            "R3",
            "alignment_base_identity",
            real_mapping["base_identity"],
            "mapping",
            "base_identity",
            "real ONT 30x through the same mapper",
        ),
        (
            "R4",
            "variant_f1",
            real_variant["f1_score"],
            "variant",
            "f1_score",
            "synthetic truth; cross-tool only and excluded from delta-to-real",
        ),
        (
            "R5",
            "assembly_reference_identity",
            real_assembly["reference_identity"],
            "assembly",
            "reference_identity",
            "real ONT 30x through the same Flye protocol",
        ),
        (
            "R5",
            "assembly_n50_bp",
            real_assembly["n50"],
            "assembly",
            "n50",
            "real ONT 30x through the same Flye protocol",
        ),
    ):
        row = {
            "layer": layer,
            "metric": metric,
            "real_anchor": real_value,
            "note": note,
        }
        for key, label in TOOLS:
            row[label] = panel[key][source_name].get(value_field, "")
        rows.append(row)

    fields = [
        "layer", "metric", "real_anchor",
        "GCerrHMM", "NanoSim", "PBSim3-sample", "badread", "note",
    ]
    csv_path = stats / "final_layer_summary_v2.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    md_path = stats / "final_layer_summary_v2.md"
    lines = [
        "# Final v2 layer summary",
        "",
        "Inputs: coverage-matched `results/framework/tables/*_30x.json`.",
        "Real anchors: 30x ONT mapping/variant tables and Flye assembly JSON.",
        "",
        "| layer | metric | real | GCerrHMM | NanoSim | PBSim3-sample | badread | note |",
        "|---|---|---:|---:|---:|---:|---:|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row['layer']} | {row['metric']} | {row['real_anchor']} | "
            f"{row['GCerrHMM']} | {row['NanoSim']} | "
            f"{row['PBSim3-sample']} | {row['badread']} | {row['note']} |"
        )
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"written: {csv_path}")
    print(f"written: {md_path}")


if __name__ == "__main__":
    main()
