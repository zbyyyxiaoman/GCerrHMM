#!/usr/bin/env python3
"""Build the six main figures from the teacher's Results-section layout."""

import argparse
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from figure_io import save_figure_pair
from project_paths import project_file


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
SPLIT_SPECIES = ["Ecoli", "Scerevisiae"]
TOOLS = [
    ("GCerrHMM", "GCerrHMM"),
    ("NanoSim", "NanoSim"),
    ("PBSim3-sample", "PBSim3\nsample"),
    ("PBSim3-errhmm", "PBSim3\nerrhmm"),
    ("badread", "badread"),
]
# Variant and assembly panels exclude PBSim3-errhmm: its FASTQ carries Q0
# placeholder quality strings, so variant calling is not evaluable, and no
# assembly run exists for it. It stays in the alignment panel, where its
# mapping statistics are valid.
TOOLS_DOWNSTREAM = [
    ("GCerrHMM", "GCerrHMM"),
    ("NanoSim", "NanoSim"),
    ("PBSim3-sample", "PBSim3\nsample"),
    ("badread", "badread"),
]
ROUTES = ["Sample-A", "GCerrHMM-1bin", "GCerrHMM"]
PALETTE = {
    "blue": "#1F4E79",
    "teal": "#2A9D8F",
    "amber": "#E9A23B",
    "red": "#D1495B",
    "purple": "#6A4C93",
    "grey": "#6B7280",
}


def read_rows(path):
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def number(value):
    try:
        if value in ("", None, "NA", "nan"):
            return math.nan
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def save_panels(fig, axes, output_dir, prefix):
    target = output_dir / "panels"
    target.mkdir(exist_ok=True)
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


def box(ax, x, y, w, h, text, color, fontsize=8):
    patch = FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.05,rounding_size=0.08",
        linewidth=0.9, edgecolor="#333333", facecolor=color,
    )
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fontsize)


def arrow(ax, start, end):
    ax.add_patch(FancyArrowPatch(
        start, end, arrowstyle="-|>", mutation_scale=11,
        linewidth=1.1, color="#444444",
    ))


def style(ax):
    ax.grid(alpha=0.16, linewidth=0.6, color="#AAB2BD")
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color("#4B5563")
        ax.spines[spine].set_linewidth(0.8)
    ax.tick_params(colors="#374151", labelsize=8)


def heatmap(ax, matrix, row_labels, col_labels, title, cmap="cividis",
            vmin=0, vmax=100, colorbar_label="Score"):
    matrix = np.asarray(matrix, dtype=float)
    if not np.isfinite(matrix).all():
        missing = [
            f"{row_labels[row]}/{col_labels[col]}"
            for row, col in zip(*np.where(~np.isfinite(matrix)))
        ]
        raise ValueError(f"{title} has missing cells: {', '.join(missing)}")
    image = ax.imshow(matrix, aspect="auto", cmap=cmap, norm=Normalize(vmin=vmin, vmax=vmax))
    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=25, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=8)
    ax.set_title(title, fontsize=10, fontweight="bold", color="#1F2937")
    ax.set_xticks(np.arange(-0.5, len(col_labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(row_labels), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.tick_params(which="major", length=3, width=0.7, colors="#374151")
    for spine in ("top", "right", "left", "bottom"):
        ax.spines[spine].set_visible(False)
    for row in range(matrix.shape[0]):
        for col in range(matrix.shape[1]):
            value = matrix[row, col]
            text = "NA" if not np.isfinite(value) else f"{value:.1f}"
            normalized = (value - vmin) / (vmax - vmin) if vmax != vmin else 0.5
            color = "#F8FAFC" if normalized > 0.62 else "#111827"
            ax.text(
                col, row, text, ha="center", va="center",
                fontsize=7.0, color=color, fontweight="normal",
            )
    return image


def metric_matrix(rows, value_field, species_list, labels):
    matrix = np.full((len(species_list), len(labels)), np.nan)
    for row in rows:
        species = row.get("species")
        label = row.get("tool_label")
        if species in species_list and label in labels:
            matrix[species_list.index(species), labels.index(label)] = number(row.get(value_field))
    return matrix


def figure1(output_dir, project_dir):
    fig, axes = plt.subplots(1, 3, figsize=(16, 5.0), constrained_layout=True)

    ax = axes[0]
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")
    box(ax, 0.2, 4.2, 1.6, 0.8, "Real ONT\nBAM", "#DCEAF7")
    box(ax, 0.2, 2.9, 1.6, 0.8, "Real HiFi\nBAM", "#DDF2D8")
    box(ax, 2.2, 3.55, 1.8, 0.9, "CIGAR error\nstate learning", "#FFF1CC")
    box(ax, 4.4, 3.55, 1.9, 0.9, "GC-conditioned\ntransition HMM", "#F7D9E3")
    box(ax, 6.7, 3.55, 1.4, 0.9, "Matched\nsimulation", "#E6DDF2")
    box(ax, 8.4, 3.55, 1.4, 0.9, "Four-layer\nevaluation", "#F4F4F4")
    arrow(ax, (1.8, 4.6), (2.2, 4.0))
    arrow(ax, (1.8, 3.3), (2.2, 3.9))
    arrow(ax, (4.0, 4.0), (4.4, 4.0))
    arrow(ax, (6.3, 4.0), (6.7, 4.0))
    arrow(ax, (8.1, 4.0), (8.4, 4.0))
    box(ax, 4.4, 1.8, 5.4, 0.8, "platform-aware evaluation: ONT | HiFi", "#FFFFFF", fontsize=8)
    ax.set_title("a  Dual-platform workflow", fontsize=10, fontweight="bold")

    ax = axes[1]
    labels = []
    values = []
    model_dir = project_dir / "data" / "trained_models"
    for species in ("Ecoli", "Scerevisiae", "Athaliana"):
        path = model_dir / f"{species}_errhmm.json"
        if not path.exists():
            continue
        import json
        model = json.loads(path.read_text(encoding="utf-8"))
        states = model["states"]
        match_index = states.index("M")
        curve = []
        for key in sorted(model["transition_probs"], key=lambda item: int(item)):
            row = model["transition_probs"][key][match_index]
            curve.append(100.0 * max(0.0, 1.0 - float(row[match_index])))
        if curve:
            labels.append(SPECIES_LABELS[species])
            values.append(curve)
    if values:
        colors = (PALETTE["blue"], PALETTE["amber"], PALETTE["teal"])
        for label, curve, color in zip(labels, values, colors):
            ax.plot(
                np.arange(len(curve)), curve, marker="o",
                linewidth=2.0, markersize=4, color=color, label=label,
            )
        ax.legend(frameon=False, fontsize=7.5)
        ax.set_xticks(np.arange(max(len(curve) for curve in values)))
        ax.set_xticklabels([f"{index}" for index in range(max(len(curve) for curve in values))])
        ax.set_xlabel("GC bin")
        ax.set_ylabel("Model-implied error (%)")
        style(ax)
    else:
        ax.axis("off")
        ax.text(0.5, 0.5, "Model matrices unavailable", ha="center", va="center")
    ax.set_title("b  GC-conditioned error state", fontsize=10, fontweight="bold")

    ax = axes[2]
    table = read_rows(
        project_file(
            project_dir,
            "docs/paper_figures/table1_species_panel.csv",
        )
    )
    if table:
        headers = ["Species", "Assembly", "Mb", "GC %", "GC SD", "ONT run"]
        rows = [
            [
                row["species"],
                row["reference_accession"],
                f"{number(row['genome_size_mb']):.1f}",
                f"{number(row['gc_percent']):.1f}",
                f"{number(row['gc_heterogeneity_sd_1kb']):.1f}",
                row["ont_read_accession"],
            ]
            for row in table
        ]
        ax.axis("off")
        table_artist = ax.table(
            cellText=rows,
            colLabels=headers,
            cellLoc="left",
            colLoc="center",
            loc="center",
            colWidths=[0.22, 0.25, 0.09, 0.09, 0.09, 0.24],
        )
        table_artist.auto_set_font_size(False)
        table_artist.set_fontsize(5.7)
        table_artist.scale(1.0, 1.35)
        for (row_index, _), cell in table_artist.get_celld().items():
            cell.set_edgecolor("#B9B9B9")
            cell.set_linewidth(0.45)
            if row_index == 0:
                cell.set_facecolor("#E8EEF5")
                cell.get_text().set_fontweight("bold")
    else:
        raise ValueError("species table unavailable")
    ax.set_title("c  Species panel and data sources", fontsize=10, fontweight="bold")

    fig.suptitle("Figure 1. GCerrHMM overview and innovation", fontsize=15, fontweight="bold")
    save_panels(fig, axes, output_dir, "figure1")
    save_figure_pair(fig, output_dir, "figure1_overview_innovation")


def figure2(output_dir, project_dir):
    reads = read_rows(project_dir / "results/framework/stats/reads_level1.csv")
    metrics = [
        ("read_length_mean", "Read length"),
        ("qv_mean", "QV"),
        ("gc_mean", "GC"),
        ("kmer_mean", "k-mer"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(16, 9.5), constrained_layout=True)
    for ax, (field, title) in zip(axes.ravel()[:4], metrics):
        matrix = metric_matrix(reads, field, EXTERNAL_SPECIES, [key for key, _ in TOOLS])
        image = heatmap(
            ax, matrix,
            [SPECIES_LABELS[item] for item in EXTERNAL_SPECIES],
            [label for _, label in TOOLS],
            title, cmap="cividis", vmin=0, vmax=100,
        )
        fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
    matrix = metric_matrix(reads, "composite_mean", EXTERNAL_SPECIES, [key for key, _ in TOOLS])
    image = heatmap(
        axes.ravel()[4], matrix,
        [SPECIES_LABELS[item] for item in EXTERNAL_SPECIES],
        [label for _, label in TOOLS],
        "Composite (ONT matched panel)", cmap="cividis", vmin=0, vmax=100,
    )
    fig.colorbar(image, ax=axes.ravel()[4], fraction=0.046, pad=0.03)
    ax = axes.ravel()[5]
    hifi = read_rows(project_dir / "results/framework/stats/hifi_level1_10x.csv")
    if hifi:
        route_labels = ["sample_A", "errhmm_1bin", "errhmm_gc"]
        display = ["Sample-A", "GCerrHMM-1bin", "GCerrHMM"]
        x = np.arange(len(route_labels))
        composite = [
            next(
                (
                    number(row.get("composite_mean", row.get("composite")))
                    for row in hifi if row.get("route") == route
                ),
                math.nan,
            )
            for route in route_labels
        ]
        kmer = [
            next(
                (
                    number(row.get("kmer_mean", row.get("kmer")))
                    for row in hifi if row.get("route") == route
                ),
                math.nan,
            )
            for route in route_labels
        ]
        composite_sd = [
            next(
                (
                    number(row.get("composite_sd", 0.0))
                    for row in hifi if row.get("route") == route
                ),
                0.0,
            )
            for route in route_labels
        ]
        kmer_sd = [
            next(
                (
                    number(row.get("kmer_sd", 0.0))
                    for row in hifi if row.get("route") == route
                ),
                0.0,
            )
            for route in route_labels
        ]
        width = 0.36
        ax.bar(
            x - width / 2, composite, width=width,
            yerr=composite_sd, capsize=3, color=PALETTE["blue"],
            label="Composite",
        )
        ax.bar(
            x + width / 2, kmer, width=width,
            yerr=kmer_sd, capsize=3, color=PALETTE["amber"], label="k-mer",
        )
        ax.set_xticks(x)
        ax.set_xticklabels(display, rotation=20, ha="right", fontsize=7.5)
        ax.set_ylim(0, 100)
        ax.set_ylabel("Score")
        ax.legend(frameon=False, fontsize=7)
        style(ax)
        ax.set_title("HiFi cross-platform check", fontsize=9.5, fontweight="bold")
    else:
        ax.axis("off")
        ax.text(0.02, 0.88, "Platform panels", ha="left", va="top", fontsize=10, fontweight="bold")
        ax.text(0.02, 0.72, "ONT: matched panel measured", ha="left", va="top", fontsize=9)
        ax.text(0.02, 0.58, "HiFi: run titan_hifi_crossplatform_one.sh", ha="left", va="top", fontsize=9)
    fig.suptitle("Figure 2. Reads-level cross-tool comparison", fontsize=15, fontweight="bold")
    save_panels(fig, axes, output_dir, "figure2")
    save_figure_pair(fig, output_dir, "figure2_reads_cross_tool")


def figure3(output_dir, project_dir):
    mapping = read_rows(project_dir / "results/framework/stats/mapping.csv")
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.2), constrained_layout=True)
    for ax, (field, title, vmin, vmax) in zip(
        axes,
        [
            ("alignment_rate_mean", "Mapping rate (%)", 80, 100),
            ("mean_mapq_mean", "Mean MAPQ", 0, 60),
            ("base_identity_mean", "Base-level identity (%)", 80, 100),
        ],
    ):
        matrix = metric_matrix(mapping, field, EXTERNAL_SPECIES, [key for key, _ in TOOLS])
        if field == "base_identity_mean":
            matrix = 100.0 * matrix
        image = heatmap(
            ax, matrix,
            [SPECIES_LABELS[item] for item in EXTERNAL_SPECIES],
            [label for _, label in TOOLS],
            title, cmap="YlGnBu", vmin=vmin, vmax=vmax,
        )
        fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
    fig.suptitle("Figure 3. Alignment / BAM-level comparison", fontsize=15, fontweight="bold")
    save_panels(fig, axes, output_dir, "figure3")
    save_figure_pair(fig, output_dir, "figure3_alignment")


def figure4(output_dir, project_dir):
    variant = read_rows(project_dir / "results/framework/stats/variant.csv")
    sv_rows = read_rows(project_dir / "results/framework/stats/sv_calling.csv")
    sv_lookup = {
        (row["species"], row["tool_label"]): row.get("f1_score", "")
        for row in sv_rows
    }
    fig, axes = plt.subplots(2, 2, figsize=(13, 9.5), constrained_layout=True)
    panels = [
        ("f1_score_mean", EXTERNAL_SPECIES, "Combined SNP + indel F1 (%)"),
        ("snp_f1_score_mean", SPLIT_SPECIES, "SNP F1 (%)"),
        ("indel_f1_score_mean", SPLIT_SPECIES, "Indel F1 (%)"),
    ]
    tool_keys = [key for key, _ in TOOLS_DOWNSTREAM]
    tool_labels = [label for _, label in TOOLS_DOWNSTREAM]
    for ax, (field, species_list, title) in zip(axes.flat, panels):
        matrix = metric_matrix(variant, field, species_list, tool_keys)
        matrix = 100.0 * matrix
        image = heatmap(
            ax, matrix,
            [SPECIES_LABELS[item] for item in species_list],
            tool_labels,
            title, cmap="cividis", vmin=0, vmax=100,
        )
        fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03)

    sv_matrix = np.full((len(SPLIT_SPECIES), len(tool_labels)), np.nan)
    for row_index, species in enumerate(SPLIT_SPECIES):
        for column_index, tool in enumerate(tool_keys):
            raw = sv_lookup.get((species, tool), "")
            if raw not in ("", None):
                sv_matrix[row_index, column_index] = 100.0 * number(raw)
    image = heatmap(
        axes[1, 1], sv_matrix,
        [SPECIES_LABELS[item] for item in SPLIT_SPECIES],
        tool_labels,
        "SV F1 (%) — spike-in", cmap="cividis", vmin=0, vmax=100,
    )
    fig.colorbar(image, ax=axes[1, 1], fraction=0.046, pad=0.03)
    fig.suptitle(
        "Figure 4. Variant-calling comparison (combined / SNP / indel / SV)",
        fontsize=15, fontweight="bold",
    )
    fig.text(
        0.5, 0.005,
        "PBSim3-errhmm is omitted (Q0 placeholder qualities make the standard "
        "caller non-evaluable); H. sapiens chr21 is shown in the combined-F1 "
        "panel because PBSim3-sample lacks a SNP/indel split.",
        ha="center", va="bottom", fontsize=8,
    )
    save_panels(fig, axes, output_dir, "figure4")
    save_figure_pair(fig, output_dir, "figure4_variant_calling")


def figure5(output_dir, project_dir):
    assembly = read_rows(project_dir / "results/framework/stats/assembly.csv")
    # Only the protocol-of-record assemblers enter the figure. Legacy raven
    # rows stay in the CSV as provenance but are not plotted as current data.
    assembly = [
        row for row in assembly
        if row.get("assembler") in ("flye", "hifiasm")
    ]
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.2), constrained_layout=True)
    for ax, (field, title, vmin, vmax) in zip(
        axes[:2],
        [
            ("n50_mean", "Assembly N50 (bp)", 0, 40_000_000),
            ("reference_identity_mean", "Assembly identity (%)", 0, 100),
        ],
    ):
        matrix = metric_matrix(assembly, field, EXTERNAL_SPECIES,
                               [key for key, _ in TOOLS_DOWNSTREAM])
        if field == "reference_identity_mean":
            matrix = 100.0 * matrix
        image = heatmap(
            ax, matrix,
            [SPECIES_LABELS[item] for item in EXTERNAL_SPECIES],
            [label for _, label in TOOLS_DOWNSTREAM],
            title, cmap="magma_r", vmin=vmin, vmax=vmax,
        )
        fig.colorbar(image, ax=ax, fraction=0.046, pad=0.03)
    axes[2].axis("off")
    phasing_results = sorted(
        project_dir.glob("results/framework/**/phasing/phasing_result.json"),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    ) or sorted(
        project_dir.glob("results/**/phasing_result.json"),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    if phasing_results:
        import json
        payload = json.loads(phasing_results[0].read_text(encoding="utf-8"))
        stats = payload.get("phase_stats", {})
        lines = [
            ("Phasing (real HiFi, whatshap)", "", True),
            ("Sample / platform",
             f"{payload.get('sample', 'NA')} / {payload.get('platform', 'NA')}", False),
            ("Mapping rate", payload.get("mapping_rate", "NA"), False),
            ("Heterozygous variants", stats.get("heterozygous_variants", "NA"), False),
            ("Phased variants", stats.get("phased", "NA"), False),
            ("Phased fraction", stats.get("phased_fraction", "NA"), False),
            ("Phase blocks", stats.get("blocks", "NA"), False),
            ("Block N50 (bp)", stats.get("block_n50", "NA"), False),
            ("Switch error", "not reported (no phased truth)", False),
        ]
        y = 0.95
        for label, value, is_header in lines:
            if is_header:
                axes[2].text(0.02, y, label, ha="left", va="top",
                             fontsize=10, fontweight="bold")
            else:
                axes[2].text(0.04, y, f"{label}:", ha="left", va="top", fontsize=8)
                axes[2].text(0.60, y, str(value), ha="left", va="top", fontsize=8)
            y -= 0.098
    else:
        axes[2].text(0.02, 0.92, "Assembly and phasing scope", ha="left", va="top", fontsize=10, fontweight="bold")
        axes[2].text(0.02, 0.76, "HiFi: hifiasm (primary)", ha="left", va="top", fontsize=9)
        axes[2].text(0.02, 0.62, "ONT: Flye (primary)", ha="left", va="top", fontsize=9)
        axes[2].text(0.02, 0.48, "Phasing: whatshap, diploid truth required", ha="left", va="top", fontsize=9)
        axes[2].text(0.02, 0.34, "Current legacy values are kept only as provenance.", ha="left", va="top", fontsize=8)
    fig.suptitle(
        "Figure 5. Assembly (Flye for ONT, hifiasm for HiFi) and phasing",
        fontsize=15, fontweight="bold",
    )
    save_panels(fig, axes, output_dir, "figure5")
    save_figure_pair(fig, output_dir, "figure5_assembly_phasing")


def figure6(output_dir, project_dir):
    stats = project_dir / "results" / "stats"
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.2), constrained_layout=True)
    gc = read_rows(stats / "profile_matched_gc_bins_seeds_20260916_010047.csv")
    order = {"1bin": 1, "5bin": 5, "10bin": 10, "20bin": 20}
    for species in ("Ecoli", "Athaliana"):
        subset = sorted([row for row in gc if row["species"] == species], key=lambda row: order[row["gc_bins"]])
        if subset:
            axes[0].errorbar(
                [order[row["gc_bins"]] for row in subset],
                [float(row["mean"]) for row in subset],
                yerr=[float(row["std"]) for row in subset],
                marker="o", capsize=4, linewidth=2.0,
                color=(PALETTE["blue"] if species == "Ecoli" else PALETTE["amber"]),
                label=SPECIES_LABELS[species],
            )
    axes[0].set_xscale("log", base=2)
    axes[0].set_xticks([1, 5, 10, 20])
    axes[0].set_xticklabels(["1", "5", "10", "20"])
    axes[0].set_xlabel("GC bins")
    axes[0].set_ylabel("Composite")
    axes[0].legend(frameon=False, fontsize=8)
    style(axes[0])
    axes[0].set_title("a  GC-bin ablation", fontsize=10, fontweight="bold")

    coverage = read_rows(stats / "coverage_replicates_20260917_142457.csv")
    axes[1].errorbar(
        np.arange(len(coverage)),
        [float(row["mean"]) for row in coverage],
        yerr=[float(row["std"]) for row in coverage],
        marker="o", capsize=4, linewidth=1.8, color="#2E8B57",
    )
    axes[1].set_xticks(np.arange(len(coverage)))
    axes[1].set_xticklabels([row["coverage"] for row in coverage])
    axes[1].set_xlabel("Training coverage")
    axes[1].set_ylabel("Composite")
    style(axes[1])
    axes[1].set_title("b  Coverage ablation", fontsize=10, fontweight="bold")

    state = read_rows(stats / "state_space_seeds_20260917_183439.csv")
    axes[2].bar(
        np.arange(len(state)),
        [float(row["composite_score_mean"]) for row in state],
        yerr=[float(row["composite_score_sd"]) for row in state],
        capsize=4, color=["#4C78A8", "#54A24B"],
    )
    axes[2].set_xticks(np.arange(len(state)))
    axes[2].set_xticklabels(["Full", "Simple"])
    axes[2].set_ylabel("Composite")
    style(axes[2])
    axes[2].set_title("c  State-space ablation", fontsize=10, fontweight="bold")
    fig.suptitle("Figure 6. Ablation and robustness", fontsize=15, fontweight="bold")
    save_panels(fig, axes, output_dir, "figure6")
    save_figure_pair(fig, output_dir, "figure6_ablation")


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
    plt.rcParams.update({
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "savefig.facecolor": "white",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "font.family": "DejaVu Sans",
        "axes.edgecolor": "#4B5563",
        "axes.linewidth": 0.8,
        "xtick.color": "#374151",
        "ytick.color": "#374151",
        "legend.frameon": False,
        "figure.dpi": 120,
    })
    figure1(output_dir, project_dir)
    figure2(output_dir, project_dir)
    figure3(output_dir, project_dir)
    figure4(output_dir, project_dir)
    figure5(output_dir, project_dir)
    figure6(output_dir, project_dir)
    print(f"GCERRHMM_MAIN_FIGURES={output_dir}")


if __name__ == "__main__":
    main()
