#!/usr/bin/env python3
"""Build background figures for the long-read simulation landscape."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

from common_io import read_csv_rows as read_csv
from figure_io import save_figure
from project_paths import project_file


PALETTE = {
    "blue": "#1F4E79",
    "teal": "#2A9D8F",
    "amber": "#E9A23B",
    "red": "#D1495B",
    "purple": "#6A4C93",
    "grey": "#6B7280",
}

def style(ax) -> None:
    ax.grid(alpha=0.16, linewidth=0.6, color="#AAB2BD")
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#4B5563")
        ax.spines[spine].set_linewidth(0.8)


def figure_b1(project_dir: Path, output_dir: Path) -> None:
    species = read_csv(
        project_file(project_dir, "docs/paper_figures/table1_species_panel.csv")
    )
    curves = read_csv(
        project_file(project_dir, "docs/paper_innovation/stats/gc_error_curve.csv")
    )
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.4), constrained_layout=True)

    ax = axes[0]
    sizes = [max(20, float(row["genome_size_mb"])) for row in species]
    colors = [
        PALETTE["blue"], PALETTE["teal"], PALETTE["amber"],
        PALETTE["red"], PALETTE["purple"], PALETTE["grey"],
    ]
    ax.scatter(
        [float(row["gc_percent"]) for row in species],
        [float(row["gc_heterogeneity_sd_1kb"]) for row in species],
        s=sizes, c=colors, alpha=0.82, edgecolor="white", linewidth=1.0,
    )
    offsets = {
        "E. coli": (0.35, 0.08),
        "S. cerevisiae": (-3.2, 0.18),
        "A. thaliana": (-7.0, -0.05),
        "D. melanogaster": (0.4, 0.08),
        "M. musculus chr19": (0.4, -0.22),
        "H. sapiens chr21": (-6.5, 0.05),
    }
    for row, color in zip(species, colors):
        name = row["species"]
        x, y = float(row["gc_percent"]), float(row["gc_heterogeneity_sd_1kb"])
        dx, dy = offsets.get(name, (0.3, 0.1))
        ax.text(x + dx, y + dy, name, fontsize=8, color=color,
                fontweight="bold")
    ax.set_xlabel("Reference GC content (%)")
    ax.set_ylabel("1-kb GC heterogeneity (SD)")
    ax.set_title("a  GC context differs across the panel", fontsize=11,
                 fontweight="bold")
    style(ax)

    ax = axes[1]
    labels = {
        "real": ("Real reads", PALETTE["blue"], "-", "o"),
        "errhmm": ("GC-aware", PALETTE["teal"], "-", "o"),
        "nanosim": ("NanoSim", PALETTE["amber"], "--", "s"),
        "pbsim": ("PBSim3", PALETTE["red"], ":", "^"),
        "badread": ("badread", PALETTE["grey"], "-.", "d"),
    }
    for dataset, (label, color, linestyle, marker) in labels.items():
        subset = sorted(
            [row for row in curves if row["dataset"] == dataset],
            key=lambda row: float(row["gc_bin"].split("-", 1)[0]),
        )
        if not subset:
            continue
        x = [float(row["gc_bin"].split("-", 1)[0]) + 5 for row in subset]
        y = [100.0 * float(row["error_rate"]) for row in subset]
        ax.plot(
            x, y, color=color, linestyle=linestyle, marker=marker,
            linewidth=2.0, markersize=4.5, label=label,
        )
    ax.set_xlabel("GC-bin midpoint (%)")
    ax.set_ylabel("Error rate (%)")
    ax.set_title("b  Example real and simulated GC-error curves",
                 fontsize=11, fontweight="bold")
    ax.legend(fontsize=8)
    style(ax)
    fig.suptitle(
        "Background Figure B1. Why local GC context is a plausible error axis",
        fontsize=14, fontweight="bold",
    )
    save_figure(fig, output_dir / "background_figure_b1_gc_context.png",
                dpi=300, bbox_inches="tight")
    save_figure(fig, output_dir / "background_figure_b1_gc_context.pdf",
                dpi=300, bbox_inches="tight")
    plt.close(fig)


def figure_b2(output_dir: Path) -> None:
    timeline = [
        (2017, "NanoSim\nread profiles", PALETTE["blue"]),
        (2018, "NPBSS\nempirical CLR", PALETTE["teal"]),
        (2022, "PBSim3\nHMM error model", PALETTE["amber"]),
        (2023, "Meta-NanoSim\nmetagenomic profile", PALETTE["purple"]),
        (2024, "Squigulator\nsignal noise", PALETTE["red"]),
        (2024, "Icarust\nadaptive sampling", PALETTE["grey"]),
        (2025, "Platinum Pedigree\nlong-read benchmark", PALETTE["blue"]),
        (2026, "Context-aware simulation\nparameter optimisation", PALETTE["teal"]),
        (2026, "NanoSimFormer\nsignal transformer", PALETTE["amber"]),
    ]
    rows = [
        "NanoSim",
        "PBSim3",
        "Squigulator",
        "NanoSimFormer",
        "Context-aware sim.",
        "GCerrHMM",
    ]
    columns = [
        "Read output",
        "Signal output",
        "Local context",
        "GC-conditioned\nerror state",
        "Multi-platform",
        "Real-data anchor",
        "Pre-registered\nevaluation",
    ]
    # 2 = primary design element; 1 = partial/supporting; 0 = not a primary
    # model axis in the cited paper/software description.
    matrix = np.array([
        [2, 0, 1, 0, 0, 0, 0],
        [2, 0, 1, 0, 2, 0, 0],
        [0, 2, 1, 0, 0, 0, 0],
        [0, 2, 2, 0, 0, 1, 0],
        [2, 0, 2, 0, 0, 2, 0],
        [2, 0, 2, 2, 2, 2, 2],
    ])
    fig, axes = plt.subplots(
        2, 1, figsize=(14, 10), constrained_layout=True,
        gridspec_kw={"height_ratios": [1, 1.25]},
    )

    ax = axes[0]
    ax.axhline(0, color="#4B5563", linewidth=1.2)
    offsets = [0.34, -0.34, 0.34, -0.34, 0.34, -0.34, 0.34, -0.34, 0.34]
    year_counts: dict[int, int] = {}
    for (year, label, color), offset in zip(timeline, offsets):
        count = year_counts.get(year, 0)
        year_counts[year] = count + 1
        jitter = (-0.16, 0.0, 0.16)[min(count, 2)]
        x_position = year + jitter
        ax.scatter(x_position, 0, s=70, color=color, zorder=3)
        ax.text(
            x_position, offset, label, ha="center",
            va="bottom" if offset > 0 else "top",
            fontsize=7.5, color="#1F2937",
        )
    ax.set_xlim(2016.5, 2026.8)
    ax.set_ylim(-0.8, 0.8)
    ax.set_yticks([])
    ax.set_xlabel("Publication year")
    ax.set_title("a  Recent long-read simulation and benchmark landscape",
                 fontsize=11, fontweight="bold")
    style(ax)

    ax = axes[1]
    display = np.full(matrix.shape, "", dtype=object)
    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            display[row, col] = {0: "–", 1: "partial", 2: "yes"}[matrix[row, col]]
    image = ax.imshow(matrix, cmap="cividis", vmin=0, vmax=2, aspect="auto")
    ax.set_xticks(np.arange(len(columns)))
    ax.set_xticklabels(columns, rotation=18, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(rows)))
    ax.set_yticklabels(rows, fontsize=8.5)
    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            value = matrix[row, col]
            color = "#F8FAFC" if value >= 1 else "#111827"
            ax.text(col, row, display[row, col], ha="center", va="center",
                    fontsize=7, color=color, fontweight="bold")
    ax.set_xticks(np.arange(-0.5, len(columns), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.9)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_title("b  Feature coverage based on primary model descriptions",
                 fontsize=11, fontweight="bold")
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.suptitle(
        "Background Figure B2. Long-read simulator landscape and evaluation gap",
        fontsize=14, fontweight="bold",
    )
    save_figure(fig, output_dir / "background_figure_b2_simulator_landscape.png",
                dpi=300, bbox_inches="tight")
    save_figure(fig, output_dir / "background_figure_b2_simulator_landscape.pdf",
                dpi=300, bbox_inches="tight")
    plt.close(fig)


def figure_b3(project_dir: Path, output_dir: Path) -> None:
    homopolymer = read_csv(
        project_file(
            project_dir,
            "docs/paper_innovation/stats/homopolymer_deletion.csv",
        )
    )
    real = sorted(
        [row for row in homopolymer if row["dataset"] == "real"],
        key=lambda row: int(row["hp_length"]),
    )
    lengths = [int(row["hp_length"]) for row in real]
    deletion_rate = [100.0 * float(row["deletion_rate"]) for row in real]

    fig, axes = plt.subplots(
        1, 2, figsize=(14, 6.2), constrained_layout=True,
        gridspec_kw={"width_ratios": [0.9, 1.1]},
    )

    ax = axes[0]
    ax.plot(
        lengths, deletion_rate, color=PALETTE["blue"], linewidth=2.2,
        marker="o", markersize=5.0,
    )
    ax.fill_between(
        lengths, deletion_rate, color=PALETTE["blue"], alpha=0.10,
    )
    ax.set_yscale("log")
    ax.set_xticks(lengths)
    ax.set_xlabel("Reference homopolymer length (bp)")
    ax.set_ylabel("Deletion rate (%)")
    ax.set_title(
        "a  Real ONT error is context-dependent",
        fontsize=11, fontweight="bold",
    )
    if len(lengths) > 1 and deletion_rate[0] > 0:
        fold_change = deletion_rate[-1] / deletion_rate[0]
        ax.annotate(
            f"{fold_change:.1f}x higher at {lengths[-1]} bp",
            xy=(lengths[-1], deletion_rate[-1]),
            xytext=(lengths[-1] - 4.0, deletion_rate[-1] * 0.38),
            fontsize=8.5, color=PALETTE["blue"], fontweight="bold",
            arrowprops={
                "arrowstyle": "-",
                "color": PALETTE["blue"],
                "linewidth": 0.9,
            },
        )
    style(ax)

    ax = axes[1]
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    layers = [
        (
            "Read-level realism",
            "length, accuracy and quality distributions",
            "NanoSim, PBSim3, 2026 benchmark [1,2,16]",
            PALETTE["blue"],
        ),
        (
            "Error-structure realism",
            "context-dependent substitutions and homopolymers",
            "2026 benchmark [16]; NanoSimFormer [15]",
            PALETTE["amber"],
        ),
        (
            "Signal-level realism",
            "raw signal and basecaller-guided generation",
            "Squigulator, Icarust, transformer models [10-12,15]",
            PALETTE["red"],
        ),
        (
            "Downstream consistency",
            "variant, assembly and metagenomic endpoints",
            "Platinum Pedigree, NanoSimFormer [13,15]",
            PALETTE["teal"],
        ),
        (
            "Instrument-claim alignment",
            "the metric must contain the property being claimed",
            "GCerrHMM, this study",
            PALETTE["purple"],
        ),
    ]
    for index, (name, detail, sources, color) in enumerate(layers):
        y = 0.86 - index * 0.18
        ax.add_patch(
            FancyBboxPatch(
                (0.02, y - 0.065), 0.96, 0.145,
                boxstyle="round,pad=0.01,rounding_size=0.02",
                facecolor=color, edgecolor=color, linewidth=1.2,
                alpha=0.13,
            )
        )
        ax.text(
            0.05, y + 0.035, name, fontsize=10, fontweight="bold",
            color=color, va="center",
        )
        ax.text(
            0.05, y - 0.015, detail, fontsize=8.5,
            color="#1F2937", va="center",
        )
        ax.text(
            0.95, y - 0.043, sources, fontsize=7.1,
            color="#4B5563", ha="right", va="center",
        )
    ax.set_title(
        "b  Recent work separates modelling and evaluation layers",
        fontsize=11, fontweight="bold",
    )

    fig.suptitle(
        "Background Figure B3. Context dependence and the dimensions of "
        "recent simulator evaluation",
        fontsize=14, fontweight="bold",
    )
    save_figure(
        fig, output_dir / "background_figure_b3_context_layers.png",
        dpi=300, bbox_inches="tight",
    )
    save_figure(
        fig, output_dir / "background_figure_b3_context_layers.pdf",
        dpi=300, bbox_inches="tight",
    )
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    project_dir = Path(args.project_dir).resolve()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    figure_b1(project_dir, output_dir)
    figure_b2(output_dir)
    figure_b3(project_dir, output_dir)
    print(f"BACKGROUND_FIGURES={output_dir}")


if __name__ == "__main__":
    main()
