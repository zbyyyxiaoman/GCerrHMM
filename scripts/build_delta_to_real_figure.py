#!/usr/bin/env python3
"""Draw the coverage-matched chr21 30x panel against the real ONT anchor.

Four panels, each with the real anchor drawn as a dashed reference line:

  (a) alignment base identity  (R3)
  (b) variant F1, overall / SNP / indel  (R4)
  (c) assembly reference identity  (R5)
  (d) assembly N50  (R5)

Every tool is simulated and evaluated at the same coverage with the same
reference, caller and assembler, so the only difference between the bars is
the simulator. The figure is descriptive only - the pre-registered ranking
that decides the headline lives in scripts/delta_to_real_decision.py.

usage: python3 build_delta_to_real_figure.py --project-dir ~/errhmm_project
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common_io import read_json
from figure_io import save_figure

SPECIES = "Hsapiens_chr21"
SUBSPECIES = "30x"
TOOLS = [
    ("errhmm", "GCerrHMM"),
    ("nanosim", "NanoSim"),
    ("pbsim", "PBSim3-sample"),
    ("badread", "badread"),
]
REAL = "real_ont30"


def find(dirs, name: str) -> dict:
    for directory in dirs:
        payload = read_json(directory / name, default={}, ignore_errors=True)
        if payload:
            return payload
    return {}


def collect(project: Path):
    dirs = [project / "results" / "framework" / "tables",
            project / "results" / "tables"]
    rows = []
    for key, label in TOOLS:
        mapping = find(dirs, f"mapping_{key}_{SPECIES}_{SUBSPECIES}.json")
        variant = find(dirs, f"variant_{key}_{SPECIES}_{SUBSPECIES}.json")
        assembly = find(dirs, f"assembly_flye_{key}_{SPECIES}_{SUBSPECIES}.json")
        if not mapping or not variant or not assembly:
            raise SystemExit(f"30x panel table missing for {label} ({key})")
        rows.append({
            "tool_key": key,
            "tool_label": label,
            "identity": mapping.get("base_identity"),
            "f1": variant.get("f1_score"),
            "snp_f1": variant.get("snp_f1_score"),
            "indel_f1": variant.get("indel_f1_score"),
            "asm_identity": assembly.get("reference_identity"),
            "asm_n50": assembly.get("n50"),
        })
    real_mapping = find(dirs, f"mapping_{REAL}_{SPECIES}.json")
    real_variant = find(dirs, f"variant_{REAL}_{SPECIES}.json")
    real_assembly = find(dirs, f"assembly_flye_{REAL}_{SPECIES}.json")
    real = {
        "identity": real_mapping.get("base_identity"),
        "f1": real_variant.get("f1_score"),
        "snp_f1": real_variant.get("snp_f1_score"),
        "indel_f1": real_variant.get("indel_f1_score"),
        "asm_identity": real_assembly.get("reference_identity"),
        "asm_n50": real_assembly.get("n50"),
    }
    missing = [k for k, v in real.items() if v is None]
    if missing:
        raise SystemExit(f"real anchor incomplete (missing {missing}); "
                         "run the anchor chain first")
    return rows, real


def bar_panel(ax, rows, field, real_value, title, ylabel, fmt="{:.3f}",
              annotate_delta=True):
    labels = [r["tool_label"] for r in rows]
    values = [r[field] if r[field] is not None else np.nan for r in rows]
    positions = np.arange(len(rows))
    colors = ["#2f6f9f" if r["tool_key"] == "errhmm" else "#b8b8b8" for r in rows]
    ax.bar(positions, values, color=colors, edgecolor="black", linewidth=0.6)
    ax.axhline(real_value, color="black", linestyle="--", linewidth=1.4,
               label="real ONT anchor")
    for pos, value in zip(positions, values):
        if np.isnan(value):
            continue
        delta = value - real_value
        text = fmt.format(value)
        if annotate_delta:
            text += f"\nΔ{delta:+.3f}"
        ax.text(pos, value, text, ha="center", va="bottom", fontsize=7)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
    ax.set_title(title, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.tick_params(axis="y", labelsize=8)
    ax.grid(axis="y", alpha=0.25, linewidth=0.5)
    ax.set_axisbelow(True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", required=True)
    parser.add_argument("--output-dir", default=None)
    args = parser.parse_args()

    project = Path(args.project_dir).expanduser()
    out_dir = Path(args.output_dir) if args.output_dir else \
        project / "docs" / "framework_figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    stats_dir = project / "results" / "framework" / "stats"
    stats_dir.mkdir(parents=True, exist_ok=True)

    rows, real = collect(project)

    fields = ["tool_key", "tool_label", "identity", "f1", "snp_f1",
              "indel_f1", "asm_identity", "asm_n50"]
    csv_path = stats_dir / "delta_to_real_panel_30x.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
        writer.writerow({"tool_key": REAL, "tool_label": "real ONT 30x",
                         **{k: real[k] for k in fields[2:]}})

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    bar_panel(axes[0][0], rows, "identity", real["identity"],
              "(a) alignment base identity", "identity")

    ax = axes[0][1]
    labels = [r["tool_label"] for r in rows]
    positions = np.arange(len(rows))
    width = 0.27
    for offset, field, colour, name in (
        (-width, "f1", "#2f6f9f", "overall"),
        (0.0, "snp_f1", "#7fb2d6", "SNP"),
        (width, "indel_f1", "#d68f7f", "indel"),
    ):
        values = [r[field] if r[field] is not None else np.nan for r in rows]
        ax.bar(positions + offset, values, width=width, color=colour,
               edgecolor="black", linewidth=0.5, label=name)
    ax.axhline(real["f1"], color="black", linestyle="--", linewidth=1.4)
    ax.axhline(real["snp_f1"], color="#7fb2d6", linestyle=":", linewidth=1.0)
    ax.axhline(real["indel_f1"], color="#d68f7f", linestyle=":", linewidth=1.0)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
    ax.set_title("(b) variant F1 (dashed = real anchor)", fontsize=10)
    ax.set_ylabel("F1", fontsize=8)
    ax.grid(axis="y", alpha=0.25, linewidth=0.5)
    ax.set_axisbelow(True)
    ax.legend(fontsize=7, frameon=False)

    bar_panel(axes[1][0], rows, "asm_identity", real["asm_identity"],
              "(c) assembly reference identity", "identity")

    ax = axes[1][1]
    values = [(r["asm_n50"] / 1e6) if r["asm_n50"] else np.nan for r in rows]
    ax.bar(positions, values,
           color=["#2f6f9f" if r["tool_key"] == "errhmm" else "#b8b8b8"
                  for r in rows],
           edgecolor="black", linewidth=0.6)
    ax.axhline(real["asm_n50"] / 1e6, color="black", linestyle="--",
               linewidth=1.4, label="real ONT anchor")
    for pos, value in zip(positions, values):
        if np.isnan(value):
            continue
        ax.text(pos, value, f"{value:.1f} Mb", ha="center", va="bottom",
                fontsize=7)
    ax.set_xticks(positions)
    ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=8)
    ax.set_title("(d) assembly N50", fontsize=10)
    ax.set_ylabel("N50 (Mb)", fontsize=8)
    ax.grid(axis="y", alpha=0.25, linewidth=0.5)
    ax.set_axisbelow(True)

    fig.suptitle(
        "Human chr21, coverage-matched 30x panel against a real ONT anchor",
        fontsize=11,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    for suffix in ("png", "pdf"):
        save_figure(
            fig,
            out_dir / f"figure_delta_to_real_30x.{suffix}",
            dpi=300,
        )
    plt.close(fig)

    print(f"written: {csv_path}")
    print(f"written: {out_dir / 'figure_delta_to_real_30x.png'}")


if __name__ == "__main__":
    main()
