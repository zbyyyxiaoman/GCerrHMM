#!/usr/bin/env python3
"""Build the teacher-framework figure suite.

Each innovation has one main multi-panel figure. The script uses only
existing frozen results plus the framework tables exported by
``export_framework_tables.py``. Missing external downstream cells are drawn
as explicit NA cells rather than silently imputed.
"""

import argparse
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from common_io import read_csv_rows as read_csv
from figure_io import save_figure_pair
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D


SPECIES = [
    "Ecoli",
    "Scerevisiae",
    "Athaliana",
    "Dmelanogaster",
    "Mmusculus_chr19",
    "Hsapiens_chr21",
]
SPECIES_LABELS = {
    "Ecoli": "E. coli",
    "Scerevisiae": "S. cerevisiae",
    "Athaliana": "A. thaliana",
    "Dmelanogaster": "D. melanogaster",
    "Mmusculus_chr19": "M. musculus chr19",
    "Hsapiens_chr21": "H. sapiens chr21",
}
EXTERNAL_SPECIES = ["Ecoli", "Scerevisiae", "Hsapiens_chr21"]
EXTERNAL_TOOLS = [
    ("errhmm_gc_ours", "GCerrHMM"),
    ("nanosim_ext", "NanoSim"),
    ("pbsim3_sample_ext", "PBSim3-sample"),
    ("badread_ext", "badread"),
]
ROUTES = [
    ("sample_A_ours", "Sample-A"),
    ("errhmm_1bin_ours", "GCerrHMM-1bin"),
    ("errhmm_gc_ours", "GCerrHMM"),
]
COLORS = {
    "Ecoli": "#4C78A8",
    "Scerevisiae": "#F58518",
    "Athaliana": "#54A24B",
    "Dmelanogaster": "#B279A2",
    "Mmusculus_chr19": "#E45756",
    "Hsapiens_chr21": "#72B7B2",
}

def as_float(value):
    try:
        if value in ("", None, "NA", "nan"):
            return math.nan
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def save_panels(fig, axes, output_dir, prefix):
    """Export every panel as a standalone PNG for editorial reuse."""
    target = panel_dir(output_dir)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    for index, axis in enumerate(np.asarray(axes, dtype=object).ravel(), start=1):
        if not axis.get_visible():
            continue
        bbox = axis.get_tightbbox(renderer).transformed(
            fig.dpi_scale_trans.inverted()
        )
        axis.figure.savefig(
            target / f"{prefix}_panel_{index:02d}.png",
            dpi=300, bbox_inches=bbox, facecolor="white",
        )


def panel_dir(output_dir):
    path = output_dir / "panels"
    path.mkdir(exist_ok=True)
    return path


def style_axis(ax):
    ax.grid(alpha=0.18, linewidth=0.6)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def annotate_heatmap(ax, matrix, fmt="{:.1f}", threshold=60):
    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            value = matrix[row, col]
            if not np.isfinite(value):
                text = "NA"
                color = "#555555"
            else:
                text = fmt.format(value)
                color = "white" if value > threshold else "#202020"
            ax.text(col, row, text, ha="center", va="center", fontsize=7.3, color=color)


def heatmap(ax, matrix, row_labels, col_labels, title, cmap="viridis", vmin=0, vmax=100):
    matrix = np.asarray(matrix, dtype=float)
    image = ax.imshow(matrix, aspect="auto", cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax))
    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=30, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=8)
    ax.set_title(title, fontsize=10.5, fontweight="bold")
    annotate_heatmap(ax, matrix)
    return image


def profile_matched_panel(ax, stats_dir):
    rows = read_csv(stats_dir / "profile_matched_gc_bins_seeds_20260916_010047.csv")
    bin_order = {"1bin": 1, "5bin": 5, "10bin": 10, "20bin": 20}
    for species in ("Ecoli", "Athaliana"):
        subset = sorted(
            [row for row in rows if row["species"] == species],
            key=lambda row: bin_order[row["gc_bins"]],
        )
        x = [bin_order[row["gc_bins"]] for row in subset]
        mean = [float(row["mean"]) for row in subset]
        sd = [float(row["std"]) for row in subset]
        ax.errorbar(
            x, mean, yerr=sd, marker="o", linewidth=2, capsize=4,
            color=COLORS[species], label=SPECIES_LABELS[species],
        )
    ax.set_xscale("log", base=2)
    ax.set_xticks([1, 5, 10, 20])
    ax.set_xticklabels(["1", "5", "10", "20"])
    ax.set_xlabel("GC bins")
    ax.set_ylabel("Level-1 composite")
    ax.set_title("Profile-matched GC-bin response", fontsize=10.5, fontweight="bold")
    ax.legend(frameon=False, fontsize=8)
    style_axis(ax)


def model_error_panel(ax, project_dir):
    model_dir = project_dir / "data" / "trained_models"
    rows = []
    labels = []
    for species in SPECIES:
        path = model_dir / f"{species}_errhmm.json"
        if not path.exists():
            continue
        payload = json.loads(path.read_text(encoding="utf-8"))
        states = payload["states"]
        match_index = states.index("M")
        values = []
        for key in sorted(payload["transition_probs"], key=lambda item: int(item)):
            row = payload["transition_probs"][key][match_index]
            values.append(100.0 * max(0.0, 1.0 - float(row[match_index])))
        rows.append(values)
        labels.append(SPECIES_LABELS[species])
    if not rows:
        ax.axis("off")
        ax.text(0.5, 0.5, "Model matrices unavailable", ha="center", va="center")
        return
    n_bins = max(len(row) for row in rows)
    matrix = np.full((len(rows), n_bins), np.nan)
    for index, row in enumerate(rows):
        matrix[index, :len(row)] = row
    image = heatmap(
        ax, matrix, labels,
        [f"{int(100 * index / n_bins):d}-{int(100 * (index + 1) / n_bins):d}%"
         for index in range(n_bins)],
        "Model-implied error rate by GC bin", cmap="magma_r", vmin=0,
        vmax=max(1.0, np.nanmax(matrix)),
    )
    return image


def route_kmer_gain_panel(ax, reads):
    internal = {
        (row["species"], row["tool_label"]): as_float(row["kmer_mean"])
        for row in reads
        if row["source"] == "internal_frozen"
    }
    values = []
    for species in SPECIES:
        b = internal.get((species, "errHMM-1bin (ours)"), math.nan)
        c = internal.get((species, "errHMM-GC (ours)"), math.nan)
        if np.isfinite(b) and np.isfinite(c):
            values.append((SPECIES_LABELS[species], c - b))
    values.sort(key=lambda item: item[1], reverse=True)
    y = np.arange(len(values))
    bars = ax.barh(
        y, [item[1] for item in values],
        color=["#2E8B57" if item[1] >= 0 else "#C44E52" for item in values],
        edgecolor="#333333", linewidth=0.5,
    )
    ax.set_yticks(y)
    ax.set_yticklabels([item[0] for item in values], fontsize=8)
    ax.axvline(0, color="#333333", linewidth=0.9)
    for bar, (_, value) in zip(bars, values):
        ax.text(
            value + (0.25 if value >= 0 else -0.25), bar.get_y() + bar.get_height() / 2,
            f"{value:+.1f}", ha="left" if value >= 0 else "right", va="center", fontsize=8,
        )
    ax.set_xlabel("k-mer subscore gain\nGC-aware minus 1-bin")
    ax.set_title("Composition gain by genome", fontsize=10.5, fontweight="bold")
    style_axis(ax)


def state_space_panel(ax, stats_dir):
    summary = read_csv(stats_dir / "state_space_seeds_20260917_183439.csv")
    detail = read_csv(stats_dir / "state_space_seeds_details_20260917_183439.csv")
    label_map = {"full_states": "Full", "simple_states": "Simple"}
    means = [float(row["composite_score_mean"]) for row in summary]
    sds = [float(row["composite_score_sd"]) for row in summary]
    x = np.arange(len(summary))
    ax.bar(x, means, yerr=sds, capsize=5, color=["#4C78A8", "#54A24B"], edgecolor="#333333", linewidth=0.5)
    for row in detail:
        index = 0 if row["configuration"] == "full_states" else 1
        ax.scatter(index, float(row["composite_score"]), color="#202020", s=24, zorder=3)
    for index, (mean, sd) in enumerate(zip(means, sds)):
        ax.text(index, mean + sd + 0.5, f"{mean:.2f}±{sd:.2f}", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([label_map[row["configuration"]] for row in summary])
    ax.set_ylabel("Level-1 composite")
    ax.set_ylim(0, max(mean + sd for mean, sd in zip(means, sds)) + 5)
    ax.set_title("State-space robustness", fontsize=10.5, fontweight="bold")
    style_axis(ax)


def coverage_panel(ax, stats_dir):
    summary = read_csv(stats_dir / "coverage_replicates_20260917_142457.csv")
    detail = read_csv(stats_dir / "coverage_replicates_details_20260917_142457.csv")
    x = np.arange(len(summary))
    means = [float(row["mean"]) for row in summary]
    sds = [float(row["std"]) for row in summary]
    ax.errorbar(x, means, yerr=sds, marker="o", linewidth=2, capsize=4, color="#2E8B57")
    coverage_index = {row["coverage"]: index for index, row in enumerate(summary)}
    for row in detail:
        ax.scatter(coverage_index[row["coverage"]], float(row["score"]), color="#202020", s=20, zorder=3)
    ax.set_xticks(x)
    ax.set_xticklabels([row["coverage"] for row in summary])
    ax.set_xlabel("Training coverage")
    ax.set_ylabel("Level-1 composite")
    ax.set_title("Coverage robustness", fontsize=10.5, fontweight="bold")
    style_axis(ax)


def improvement_panel(ax, project_dir):
    candidates = []
    for pattern in (
        "results/gc_improvement*/**/gc_improvement_summary.json",
        "results/gc_improvement*/**/gc_improvement_summary.json",
    ):
        candidates.extend(project_dir.glob(pattern))
    candidates = sorted(set(candidates), key=lambda path: path.stat().st_mtime, reverse=True)
    if not candidates:
        ax.axis("off")
        ax.text(
            0.5, 0.5,
            "Run reproduce.sh --gc-demo\nto generate the C-minus-B gate",
            ha="center", va="center", fontsize=10,
        )
        ax.set_title("Reproducible GC-aware gain", fontsize=10.5, fontweight="bold")
        return
    payload = json.loads(candidates[0].read_text(encoding="utf-8"))
    deltas = payload["gc_minus_1bin"]
    metrics = ["composite", "kmer", "gc", "read_length", "qv"]
    labels = ["Composite", "k-mer", "GC", "Read length", "QV"]
    values = [float(deltas.get(metric, 0.0)) for metric in metrics]
    colors = ["#2E8B57" if value >= 0 else "#C44E52" for value in values]
    x = np.arange(len(values))
    ax.bar(x, values, color=colors, edgecolor="#333333", linewidth=0.5)
    ax.axhline(0, color="#222222", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("GC-aware minus 1-bin")
    ax.set_title("Reproducible GC-aware gain", fontsize=10.5, fontweight="bold")
    for index, value in enumerate(values):
        offset = 0.08 if value >= 0 else -0.08
        ax.text(index, value + offset, f"{value:+.2f}", ha="center", va="bottom" if value >= 0 else "top", fontsize=8)
    style_axis(ax)


def build_i1(output_dir, project_dir):
    stats_dir = project_dir / "results" / "stats"
    reads = read_csv(project_dir / "results/framework/stats/reads_level1.csv")
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 9.2), constrained_layout=True)
    profile_matched_panel(axes[0, 0], stats_dir)
    route_kmer_gain_panel(axes[0, 1], reads)
    heat_image = model_error_panel(axes[0, 2], project_dir)
    if heat_image is not None:
        fig.colorbar(heat_image, ax=axes[0, 2], fraction=0.046, pad=0.03, label="Error (%)")
    state_space_panel(axes[1, 0], stats_dir)
    coverage_panel(axes[1, 1], stats_dir)
    improvement_panel(axes[1, 2], project_dir)

    fig.suptitle(
        "Innovation I: GC-conditioned error modeling",
        fontsize=17, fontweight="bold",
    )
    save_panels(fig, axes, output_dir, "figure_i1")
    save_figure_pair(fig, output_dir, "figure_i1_gc_aware_model")


def external_matrix(rows, value_field, value_key, species_list, tools):
    matrix = np.full((len(species_list), len(tools)), np.nan)
    for row in rows:
        if row["species"] not in species_list or row["tool_label"] not in [label for _, label in tools]:
            continue
        species_index = species_list.index(row["species"])
        tool_index = [label for _, label in tools].index(row["tool_label"])
        matrix[species_index, tool_index] = as_float(row[value_field])
    return matrix


def build_i2(output_dir, project_dir):
    reads = read_csv(project_dir / "results/framework/stats/reads_level1.csv")
    tools = EXTERNAL_TOOLS
    fig, axes = plt.subplots(2, 3, figsize=(16, 9.6), constrained_layout=True)

    composite = external_matrix(reads, "composite_mean", "composite", EXTERNAL_SPECIES, tools)
    image = heatmap(
        axes[0, 0], composite,
        [SPECIES_LABELS[item] for item in EXTERNAL_SPECIES],
        [label for _, label in tools],
        "Cross-tool composite", cmap="YlGnBu", vmin=0, vmax=100,
    )
    fig.colorbar(image, ax=axes[0, 0], fraction=0.046, pad=0.03)

    metric_names = ["read_length", "qv", "gc", "kmer"]
    metric_labels = ["Read length", "QV", "GC", "k-mer"]
    metric_matrix = []
    metric_rows = []
    for species in EXTERNAL_SPECIES:
        for _, label in tools:
            row = next(
                (item for item in reads if item["species"] == species and item["tool_label"] == label),
                None,
            )
            metric_rows.append(f"{SPECIES_LABELS[species]} | {label}")
            metric_matrix.append([
                as_float(row[f"{metric}_mean"]) if row else math.nan
                for metric in metric_names
            ])
    image = heatmap(
        axes[0, 1], metric_matrix, metric_rows, metric_labels,
        "Read-level subscores", cmap="viridis", vmin=0, vmax=100,
    )
    fig.colorbar(image, ax=axes[0, 1], fraction=0.046, pad=0.03)

    ranking_matrix = np.full_like(composite, np.nan)
    for row_index in range(composite.shape[0]):
        order = np.argsort(-composite[row_index])
        for rank, col_index in enumerate(order, start=1):
            ranking_matrix[row_index, col_index] = rank
    heatmap(
        axes[0, 2], ranking_matrix,
        [SPECIES_LABELS[item] for item in EXTERNAL_SPECIES],
        [label for _, label in tools],
        "Within-species rank (1 = best)", cmap="RdYlGn_r", vmin=1, vmax=4,
    )

    pca_path = project_dir / "docs/paper_innovation/stats/kmer_pca.csv"
    pca = read_csv(pca_path) if pca_path.exists() else []
    ax = axes[1, 0]
    if pca:
        labels = ["Real", "Route A", "Route B", "Route C", "errHMM", "NanoSim", "PBSim3"]
        palette = plt.get_cmap("tab10")
        for index, label in enumerate(labels):
            subset = [row for row in pca if row["dataset"] == label]
            if not subset:
                continue
            ax.scatter(
                [float(row["pc1"]) for row in subset],
                [float(row["pc2"]) for row in subset],
                s=24, alpha=0.75, color=palette(index), label=label,
            )
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.set_title("k-mer composition PCA", fontsize=10.5, fontweight="bold")
        ax.legend(fontsize=7, frameon=False, ncol=2)
        style_axis(ax)
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "PCA data unavailable", ha="center", va="center")

    ax = axes[1, 1]
    species_x = np.arange(len(EXTERNAL_SPECIES))
    width = 0.18
    for tool_index, (tool_key, label) in enumerate(tools):
        values = [
            next(
                (as_float(row["composite_mean"]) for row in reads
                 if row["species"] == species and row["tool_label"] == label),
                math.nan,
            )
            for species in EXTERNAL_SPECIES
        ]
        ax.bar(
            species_x + (tool_index - 1.5) * width, values, width=width,
            label=label, color=plt.get_cmap("tab10")(tool_index),
        )
    ax.set_xticks(species_x)
    ax.set_xticklabels([SPECIES_LABELS[item] for item in EXTERNAL_SPECIES], rotation=20, ha="right")
    ax.set_ylabel("Composite")
    ax.set_title("Composite is species-dependent", fontsize=10.5, fontweight="bold")
    ax.legend(fontsize=7, frameon=False, ncol=2)
    style_axis(ax)

    status = read_csv(project_dir / "results/framework/stats/four_layer_status.csv")
    layer_keys = ["reads", "mapping", "variant", "assembly"]
    selected = [
        row for row in status
        if row["species"] in EXTERNAL_SPECIES and row["tool_label"] in [label for _, label in tools]
    ]
    status_matrix = np.asarray([
        [1.0 if row[layer] == "measured" else 0.0 for layer in layer_keys]
        for row in selected
    ])
    image = heatmap(
        axes[1, 2], status_matrix,
        [f"{SPECIES_LABELS[row['species']]} | {row['tool_label']}" for row in selected],
        ["Reads", "Mapping", "Variant", "Assembly"],
        "Evaluation coverage", cmap="RdYlGn", vmin=0, vmax=1,
    )
    for row in range(status_matrix.shape[0]):
        for col in range(status_matrix.shape[1]):
            axes[1, 2].text(
                col, row, "measured" if status_matrix[row, col] else "NA",
                ha="center", va="center", fontsize=6.5,
                color="#202020",
            )

    fig.suptitle(
        "Innovation II: matched-data cross-tool benchmarking",
        fontsize=17, fontweight="bold",
    )
    save_panels(fig, axes, output_dir, "figure_i2")
    save_figure_pair(fig, output_dir, "figure_i2_cross_tool_benchmark")


def layer_heatmap_panel(ax, rows, value_field, title, cmap, vmin=0, vmax=100):
    labels = []
    matrix = []
    for species in EXTERNAL_SPECIES:
        for tool_key, tool_label in EXTERNAL_TOOLS:
            match = next(
                (row for row in rows if row["species"] == species and row["tool_label"] == tool_label),
                None,
            )
            labels.append(f"{SPECIES_LABELS[species]} | {tool_label}")
            matrix.append([as_float(match[value_field]) if match else math.nan])
    image = heatmap(ax, matrix, labels, [title], title, cmap=cmap, vmin=vmin, vmax=vmax)
    return image


def build_i3(output_dir, project_dir):
    mapping = read_csv(project_dir / "results/framework/stats/mapping.csv")
    variant = read_csv(project_dir / "results/framework/stats/variant.csv")
    assembly = read_csv(project_dir / "results/framework/stats/assembly.csv")
    reads = read_csv(project_dir / "results/framework/stats/reads_level1.csv")
    status = read_csv(project_dir / "results/framework/stats/four_layer_status.csv")

    fig, axes = plt.subplots(2, 3, figsize=(16.4, 10), constrained_layout=True)

    map_matrix = np.full((len(EXTERNAL_SPECIES), len(EXTERNAL_TOOLS)), np.nan)
    for row in mapping:
        if row["species"] in EXTERNAL_SPECIES and row["tool_label"] in [x[1] for x in EXTERNAL_TOOLS]:
            map_matrix[
                EXTERNAL_SPECIES.index(row["species"]),
                [x[1] for x in EXTERNAL_TOOLS].index(row["tool_label"]),
            ] = as_float(row["alignment_rate_mean"])
    image = heatmap(
        axes[0, 0], map_matrix,
        [SPECIES_LABELS[item] for item in EXTERNAL_SPECIES],
        [label for _, label in EXTERNAL_TOOLS],
        "Mapping rate (%)", cmap="YlGnBu", vmin=0, vmax=100,
    )
    fig.colorbar(image, ax=axes[0, 0], fraction=0.046, pad=0.03)

    variant_matrix = np.full((len(EXTERNAL_SPECIES), len(EXTERNAL_TOOLS)), np.nan)
    for row in variant:
        if row["species"] in EXTERNAL_SPECIES and row["tool_label"] in [x[1] for x in EXTERNAL_TOOLS]:
            variant_matrix[
                EXTERNAL_SPECIES.index(row["species"]),
                [x[1] for x in EXTERNAL_TOOLS].index(row["tool_label"]),
            ] = 100.0 * as_float(row["f1_score_mean"])
    image = heatmap(
        axes[0, 1], variant_matrix,
        [SPECIES_LABELS[item] for item in EXTERNAL_SPECIES],
        [label for _, label in EXTERNAL_TOOLS],
        "Variant F1 (%)", cmap="viridis", vmin=0, vmax=100,
    )
    fig.colorbar(image, ax=axes[0, 1], fraction=0.046, pad=0.03)

    assembly_matrix = np.full((len(EXTERNAL_SPECIES), len(EXTERNAL_TOOLS)), np.nan)
    for row in assembly:
        if row["species"] in EXTERNAL_SPECIES and row["tool_label"] in [x[1] for x in EXTERNAL_TOOLS]:
            assembly_matrix[
                EXTERNAL_SPECIES.index(row["species"]),
                [x[1] for x in EXTERNAL_TOOLS].index(row["tool_label"]),
            ] = 100.0 * as_float(row["reference_identity_mean"])
    image = heatmap(
        axes[0, 2], assembly_matrix,
        [SPECIES_LABELS[item] for item in EXTERNAL_SPECIES],
        [label for _, label in EXTERNAL_TOOLS],
        "Assembly identity (%)", cmap="magma_r", vmin=0, vmax=100,
    )
    fig.colorbar(image, ax=axes[0, 2], fraction=0.046, pad=0.03)

    ax = axes[1, 0]
    internal_mapping = [
        row for row in mapping
        if row["source"] == "internal_frozen" and row["species"] in EXTERNAL_SPECIES
    ]
    internal_variant = [
        row for row in variant
        if row["source"] == "internal_frozen" and row["species"] in EXTERNAL_SPECIES
    ]
    internal_assembly = [
        row for row in assembly
        if row["source"] == "internal_frozen" and row["species"] in EXTERNAL_SPECIES
    ]
    route_labels = [label for _, label in ROUTES]
    route_matrix = np.full((len(EXTERNAL_SPECIES), len(route_labels)), np.nan)
    for row in internal_variant:
        if row["tool_label"] in route_labels:
            route_matrix[
                EXTERNAL_SPECIES.index(row["species"]),
                route_labels.index(row["tool_label"]),
            ] = 100.0 * as_float(row["f1_score_mean"])
    heatmap(
        ax, route_matrix,
        [SPECIES_LABELS[item] for item in EXTERNAL_SPECIES],
        route_labels,
        "Internal-route variant F1 (%)", cmap="viridis", vmin=0, vmax=100,
    )

    ax = axes[1, 1]
    points = []
    for row in reads:
        if row["species"] not in EXTERNAL_SPECIES:
            continue
        f1_row = next(
            (item for item in variant
             if item["species"] == row["species"] and item["tool_label"] == row["tool_label"]),
            None,
        )
        if f1_row is None:
            continue
        points.append((
            as_float(row["composite_mean"]),
            100.0 * as_float(f1_row["f1_score_mean"]),
            row["tool_label"],
            row["species"],
        ))
    for x, y, tool_label, species in points:
        ax.scatter(
            x, y, s=42, color=COLORS.get(species, "#555555"),
            edgecolor="#202020", linewidth=0.4, alpha=0.85,
        )
        ax.annotate(tool_label.split(" (")[0], (x, y), xytext=(3, 3), textcoords="offset points", fontsize=6.5)
    if len(points) >= 2:
        x = [point[0] for point in points]
        y = [point[1] for point in points]
        if len(set(x)) > 1:
            correlation = float(np.corrcoef(x, y)[0, 1])
            ax.text(0.03, 0.96, f"Pearson r = {correlation:.2f}", transform=ax.transAxes, va="top", fontsize=8)
    ax.set_xlabel("Read-level composite")
    ax.set_ylabel("Variant F1 (%)")
    ax.set_title("Read fidelity vs downstream utility", fontsize=10.5, fontweight="bold")
    style_axis(ax)

    ax = axes[1, 2]
    ax.axis("off")
    layer_rows = [row for row in status if row["species"] in EXTERNAL_SPECIES]
    measured = {
        layer: sum(row[layer] == "measured" for row in layer_rows)
        for layer in ("reads", "mapping", "variant", "assembly")
    }
    total = len(layer_rows)
    lines = ["Four-layer coverage", ""]
    for layer, label in (
        ("reads", "1. Reads / FASTQ"),
        ("mapping", "2. Alignment / BAM"),
        ("variant", "3. Variant calling"),
        ("assembly", "4. Assembly"),
    ):
        lines.append(f"{label}: {measured[layer]}/{total} measured")
    lines += ["", "Phasing is not claimed here.", "HiFi sources are validated separately."]
    for index, line in enumerate(lines):
        ax.text(
            0.02, 0.94 - index * 0.11, line, ha="left", va="top",
            fontsize=11 if index == 0 else 9,
            fontweight="bold" if index == 0 else "normal",
        )
    ax.set_title("Evidence coverage", fontsize=10.5, fontweight="bold")

    fig.suptitle(
        "Innovation III: task-aware evaluation beyond read distributions",
        fontsize=17, fontweight="bold",
    )
    save_panels(fig, axes, output_dir, "figure_i3")
    save_figure_pair(fig, output_dir, "figure_i3_task_aware_evaluation")


def write_captions(output_dir):
    text = """# Teacher-framework figure captions

## Figure I1. GC-conditioned error modeling

Main evidence for innovation I. Panel A shows the profile-matched
GC-bin sweep for *E. coli* and *A. thaliana* (mean ± SD, three seeds).
Panel B shows the species-wise k-mer subscore gain of the GC-aware model
over its one-bin control. Panel C summarizes the error probability implied
by each trained transition matrix across GC bins. Panels D and E show
state-space and training-coverage robustness. The lower-right panel records
the provenance boundary between the internal six-species panel and the
matched external panel.

## Figure I2. Matched-data cross-tool benchmarking

Main evidence for innovation II. Panels A-C compare composite scores,
read-level subscores, and within-species rankings for errHMM-GC (ours),
NanoSim, PBSim3-sample, and badread on the three-species matched panel.
Panel D is the k-mer composition PCA. Panel E shows that composite
rankings are species-dependent. Panel F explicitly marks which downstream
layers have been measured for each species and tool; NA means not yet
measured, not zero.

## Figure I3. Task-aware evaluation beyond read distributions

Main evidence for innovation III. The upper row maps the external
three-species panel through alignment, variant calling, and assembly.
The lower-left panel summarizes the internal three-route variant response.
The lower-middle panel tests whether read-level fidelity tracks downstream
variant utility. The lower-right panel provides a complete measured/NA
coverage audit. Phasing is deliberately excluded until a diploid,
truth-annotated HiFi dataset is run through the same pipeline.
"""
    (output_dir / "figure_captions.md").write_text(text, encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    project_dir = Path(args.project_dir).resolve()
    output_dir = Path(args.output_dir)
    if output_dir.exists() and any(output_dir.iterdir()):
        raise SystemExit(f"refusing to overwrite existing output directory: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    panel_dir(output_dir)
    plt.rcParams.update({
        "font.size": 9,
        "axes.titlesize": 10.5,
        "axes.labelsize": 9,
        "savefig.facecolor": "white",
    })
    build_i1(output_dir, project_dir)
    build_i2(output_dir, project_dir)
    build_i3(output_dir, project_dir)
    write_captions(output_dir)
    print(f"TEACHER_FRAMEWORK_FIGURES={output_dir}")


if __name__ == "__main__":
    main()
