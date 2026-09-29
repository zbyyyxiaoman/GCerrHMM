#!/usr/bin/env python3
"""Export the teacher-facing four-layer evaluation matrix."""

import argparse
import csv
import glob
import json
import math
import statistics
from pathlib import Path


SPECIES = (
    "Ecoli",
    "Scerevisiae",
    "Athaliana",
    "Dmelanogaster",
    "Mmusculus_chr19",
    "Hsapiens_chr21",
)
EXTERNAL_SPECIES = ("Ecoli", "Scerevisiae", "Hsapiens_chr21")
INTERNAL_TOOLS = {
    "route_A_sample": "Sample-A",
    "route_B_qshmm": "GCerrHMM-1bin",
    "route_C_errhmm": "GCerrHMM",
}
EXTERNAL_TOOLS = {
    "errhmm": "GCerrHMM",
    "nanosim": "NanoSim",
    "pbsim": "PBSim3-sample",
    "pbsim_errhmm": "PBSim3-errhmm",
    "badread": "badread",
}

# Assemblers preferred for the frozen external protocol, in priority order.
ASSEMBLER_PREFIXES = ("flye", "hifiasm")


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def table_dirs(project):
    """Result-table roots in priority order.

    The teacher-framework runs write into ``results/framework/tables``; the
    earlier phase-5 pipeline wrote into ``results/tables``. Both are searched
    so one export covers the whole project.
    """
    return [
        project / "results" / "framework" / "tables",
        project / "results" / "tables",
    ]


def find_table(dirs, name):
    for directory in dirs:
        candidate = directory / name
        if candidate.exists() and candidate.stat().st_size > 0:
            return candidate
    return None


def assembler_from_source(source_name: str) -> str:
    """Recover which assembler produced a row from its file name.

    Rows from `internal_frozen` predate the platform-matched protocol and were
    produced with raven, so they are labelled as legacy rather than silently
    mixed with Flye/hifiasm results.
    """
    name = (source_name or "").lower()
    if name.startswith("assembly_flye"):
        return "flye"
    if name.startswith("assembly_hifiasm"):
        return "hifiasm"
    if name == "internal_frozen" or name.startswith("assembly_"):
        return "raven (legacy)"
    return "unknown"


def real_species_from_source(stem: str) -> str:
    """Recover the species from a real-data assembly file name.

    Examples:
      assembly_flye_real_ont_Hsapiens_chr21_fair_v2 -> Hsapiens_chr21
      assembly_hifiasm_hifi_Hsapiens_chr21_fair_v2  -> Hsapiens_chr21
    """
    body = stem.replace("assembly_", "", 1)
    for marker in ("real_ont_", "hifi_"):
        if marker in body:
            body = body.split(marker, 1)[1]
            break
    body = body.rsplit("_fair_v2", 1)[0]
    return body or "unknown"


def mean_sd(values):
    values = [float(value) for value in values if value not in ("", None)]
    if not values:
        return "", ""
    if len(values) == 1:
        return values[0], 0.0
    return statistics.fmean(values), statistics.stdev(values)


def wilson_ci(successes: int, total: int, z: float = 1.959963985) -> str:
    """Wilson score interval as a compact 'low-high' string."""
    if total <= 0:
        return ""
    phat = successes / total
    denom = 1 + z * z / total
    centre = phat + z * z / (2 * total)
    margin = z * math.sqrt(phat * (1 - phat) / total + z * z / (4 * total * total))
    low = max(0.0, (centre - margin) / denom)
    high = min(1.0, (centre + margin) / denom)
    return f"{low:.4f}-{high:.4f}"


def internal_rows(dirs, prefix):
    rows = {}
    for species in SPECIES:
        for route, label in INTERNAL_TOOLS.items():
            paths = []
            for directory in dirs:
                paths = sorted(directory.glob(f"{prefix}_{species}_{route}_r*.json"))
                if paths:
                    break
            if not paths:
                continue
            payloads = [read_json(path) for path in paths]
            rows[(species, label)] = {
                "species": species,
                "tool_key": route,
                "tool_label": label,
                "source": "internal_frozen",
                "n": len(payloads),
                "payloads": payloads,
            }
    return rows


def external_rows(dirs, prefix):
    rows = {}
    for species in EXTERNAL_SPECIES:
        for tool, label in EXTERNAL_TOOLS.items():
            # Preferred suffixes first: the fair/alignment re-runs supersede
            # the phase-5 originals, which are kept only as a fallback.
            if species == "Hsapiens_chr21" and prefix == "mapping":
                suffixes = (
                    "_30x.json", "_alignment.json",
                    "_fair_v2.json", "_fair.json",
                )
            elif species == "Hsapiens_chr21" and prefix == "variant":
                suffixes = (
                    "_30x.json", "_fair_v3.json",
                    "_fair_v2.json", "_fair.json",
                )
            else:
                suffixes = (
                    ("_alignment.json", "_fair_v2.json", "_fair.json")
                    if prefix == "mapping"
                    else ("_fair_v3.json", "_fair_v2.json", "_fair.json")
                )
            names = [f"{prefix}_{tool}_{species}{suffix}" for suffix in suffixes]
            names.append(f"{prefix}_{tool}_{species}.json")
            if prefix == "assembly":
                if species == "Hsapiens_chr21":
                    names = [
                        f"assembly_flye_{tool}_{species}_30x.json",
                        f"assembly_hifiasm_{tool}_{species}_30x.json",
                    ] + [
                        f"assembly_{assembler}_{tool}_{species}_fair_v2.json"
                        for assembler in ASSEMBLER_PREFIXES
                    ] + [f"assembly_{tool}_{species}_fair_v2.json"]
                else:
                    names = [
                        f"assembly_{assembler}_{tool}_{species}_fair_v2.json"
                        for assembler in ASSEMBLER_PREFIXES
                    ] + [f"assembly_{tool}_{species}_fair_v2.json"]
            if prefix == "cross_tool_fixed" and tool == "pbsim":
                names.append(f"cross_tool_fixed_pbsim_{species}_full.json")
            path = next(
                (hit for name in names if (hit := find_table(dirs, name))), None
            )
            if (
                species == "Hsapiens_chr21"
                and prefix in {"mapping", "variant", "assembly"}
                and tool != "pbsim_errhmm"
                and (path is None or "_30x.json" not in path.name)
            ):
                raise SystemExit(
                    f"{prefix}: Hsapiens_chr21 must use the coverage-matched "
                    f"30x file; got {None if path is None else path.name}"
                )
            if path is None:
                continue
            rows[(species, label)] = {
                "species": species,
                "tool_key": tool,
                "tool_label": label,
                "source": path.name,
                "n": 1,
                "payloads": [read_json(path)],
            }
    return rows


def write_csv(path, fieldnames, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def build_reads_table(project):
    dirs = table_dirs(project)
    route_rows = internal_rows(dirs, "level1")
    external_rows_by_key = external_rows(dirs, "cross_tool_fixed")
    rows = []
    for key in sorted(set(route_rows) | set(external_rows_by_key), key=str):
        source = route_rows.get(key) or external_rows_by_key[key]
        payloads = source["payloads"]
        composites = []
        metrics = {name: [] for name in ("read_length", "qv", "gc", "kmer")}
        for payload in payloads:
            composite = payload.get("composite_score", {})
            composites.append(float(composite.get("composite_score", 0.0)))
            subs = composite.get("sub_scores", {})
            for name in metrics:
                if name in subs:
                    metrics[name].append(subs[name])
        composite_mean, composite_sd = mean_sd(composites)
        row = {
            "species": source["species"],
            "tool_key": source["tool_key"],
            "tool_label": source["tool_label"],
            "source": source["source"],
            "n": source["n"],
            "composite_mean": composite_mean,
            "composite_sd": composite_sd,
        }
        for name, values in metrics.items():
            mean, sd = mean_sd(values)
            row[f"{name}_mean"] = mean
            row[f"{name}_sd"] = sd
        rows.append(row)
    fields = [
        "species", "tool_key", "tool_label", "source", "n",
        "composite_mean", "composite_sd",
        "read_length_mean", "read_length_sd",
        "qv_mean", "qv_sd",
        "gc_mean", "gc_sd",
        "kmer_mean", "kmer_sd",
    ]
    write_csv(project / "results/framework/stats/reads_level1.csv", fields, rows)
    return rows


def build_mapping_table(project):
    dirs = table_dirs(project)
    internal = internal_rows(dirs, "mapping")
    external = external_rows(dirs, "mapping")
    rows = []
    alignment_metrics = (
        "base_identity", "mismatch_rate", "insertion_rate", "deletion_rate",
        "mean_read_identity", "median_read_identity",
        "reads_identity_ge_90", "reads_identity_ge_95",
        "mean_aligned_length",
    )
    for key in sorted(set(internal) | set(external), key=str):
        # The mapping layer is reported from the fair cross-tool panel when
        # available, so every tool is measured on the same simulation protocol.
        source = external.get(key) or internal.get(key) or external[key]
        payloads = source["payloads"]
        alignment_rate = mean_sd([
            payload.get("alignment_rate") for payload in payloads
            if payload.get("alignment_rate") is not None
        ])
        mean_mapq = mean_sd([
            payload.get("mean_mapq") for payload in payloads
            if payload.get("mean_mapq") is not None
        ])
        total_reads = mean_sd([
            payload.get("total_reads", 0) for payload in payloads
        ])
        row = {
            "species": source["species"],
            "tool_key": source["tool_key"],
            "tool_label": source["tool_label"],
            "source": source["source"],
            "n": source["n"],
            "alignment_rate_mean": alignment_rate[0],
            "alignment_rate_sd": alignment_rate[1],
            "mean_mapq_mean": mean_mapq[0],
            "mean_mapq_sd": mean_mapq[1],
            "total_reads_mean": total_reads[0],
        }
        for metric in alignment_metrics:
            mean, sd = mean_sd([
                payload.get(metric) for payload in payloads
                if payload.get(metric) is not None
            ])
            row[f"{metric}_mean"] = mean
            row[f"{metric}_sd"] = sd
        rows.append(row)
    fields = [
        "species", "tool_key", "tool_label", "source", "n",
        "alignment_rate_mean", "alignment_rate_sd",
        "mean_mapq_mean", "mean_mapq_sd", "total_reads_mean",
    ]
    for metric in alignment_metrics:
        fields += [f"{metric}_mean", f"{metric}_sd"]
    write_csv(project / "results/framework/stats/mapping.csv", fields, rows)
    return rows


def build_variant_table(project):
    dirs = table_dirs(project)
    internal = internal_rows(dirs, "variant")
    external = external_rows(dirs, "variant")
    rows = []
    for key in sorted(set(internal) | set(external), key=str):
        source = external.get(key) or internal.get(key) or external[key]
        row = {
            "species": source["species"],
            "tool_key": source["tool_key"],
            "tool_label": source["tool_label"],
            "source": source["source"],
            "n": source["n"],
        }
        # PBSim3's native errhmm mode writes '!' (Phred 0) for every base, so a
        # quality-aware caller discards the whole pileup and returns an empty
        # VCF. Reporting 0.0 would read as "calls nothing"; it is not evaluable
        # under a quality filter, and that is recorded explicitly instead.
        q0_placeholder = source["tool_key"] == "pbsim_errhmm" and all(
            payload.get("total_variants") == 0 for payload in source["payloads"]
        )
        row["note"] = (
            "not evaluable: input FASTQ carries Q0 placeholder qualities"
            if q0_placeholder
            else ""
        )
        for metric in (
            "precision", "recall", "f1_score",
            "snp_precision", "snp_recall", "snp_f1_score",
            "indel_precision", "indel_recall", "indel_f1_score",
            "snp_count", "indel_count", "total_variants",
        ):
            if q0_placeholder:
                row[f"{metric}_mean"] = "NA"
                row[f"{metric}_sd"] = ""
                continue
            values = [payload.get(metric) for payload in source["payloads"] if payload.get(metric) is not None]
            if metric in ("snp_count", "indel_count", "total_variants"):
                mean, sd = mean_sd(values)
                row[f"{metric}_mean"] = mean
                row[f"{metric}_sd"] = sd
            else:
                mean, sd = mean_sd(values)
                row[f"{metric}_mean"] = mean
                row[f"{metric}_sd"] = sd
        rows.append(row)
    fields = ["species", "tool_key", "tool_label", "source", "n", "note"]
    for metric in (
        "precision", "recall", "f1_score",
        "snp_precision", "snp_recall", "snp_f1_score",
        "indel_precision", "indel_recall", "indel_f1_score",
        "snp_count", "indel_count", "total_variants",
    ):
        fields += [f"{metric}_mean", f"{metric}_sd"]
    write_csv(project / "results/framework/stats/variant.csv", fields, rows)
    return rows


def build_assembly_table(project):
    dirs = table_dirs(project)
    internal = internal_rows(dirs, "assembly")
    external = external_rows(dirs, "assembly")
    rows = []
    for key in sorted(set(internal) | set(external), key=str):
        source = external.get(key) or internal.get(key) or external[key]
        row = {
            "species": source["species"],
            "tool_key": source["tool_key"],
            "tool_label": source["tool_label"],
            "source": source["source"],
            "n": source["n"],
        }
        row["assembler"] = assembler_from_source(source["source"])
        for metric in ("num_contigs", "total_length", "max_contig", "n50", "mean_contig", "reference_identity", "total_aligned_bases"):
            values = [payload.get(metric) for payload in source["payloads"] if payload.get(metric) is not None]
            mean, sd = mean_sd(values)
            row[f"{metric}_mean"] = mean
            row[f"{metric}_sd"] = sd
        rows.append(row)
    # Real-data reference points. These are the controls the simulated
    # assemblies must be compared against on the *same* platform and assembler:
    # real ONT reads through Flye, and real HiFi reads through hifiasm.
    for pattern, key, label, assembler in (
        ("assembly_flye_real_ont_*_fair_v2.json", "ont_real", "Real ONT", "flye"),
        ("assembly_hifiasm_hifi_*_fair_v2.json", "hifi_real", "Real HiFi", "hifiasm"),
    ):
        paths = sorted(
            p for directory in dirs for p in directory.glob(pattern)
        )
        for path in paths:
            payload = read_json(path)
            row = {
                "species": real_species_from_source(path.stem),
                "tool_key": key,
                "tool_label": label,
                "source": path.name,
                "assembler": assembler,
                "n": 1,
            }
            for metric in (
                "num_contigs", "total_length", "max_contig", "n50",
                "mean_contig", "reference_identity", "total_aligned_bases",
            ):
                value = payload.get(metric)
                row[f"{metric}_mean"] = value if value is not None else ""
                row[f"{metric}_sd"] = 0.0
            rows.append(row)
    fields = ["species", "tool_key", "tool_label", "assembler", "source", "n"]
    for metric in ("num_contigs", "total_length", "max_contig", "n50", "mean_contig", "reference_identity", "total_aligned_bases"):
        fields += [f"{metric}_mean", f"{metric}_sd"]
    write_csv(project / "results/framework/stats/assembly.csv", fields, rows)
    return rows


def build_coverage_matched_table(project):
    """Coverage-matched chr21 assembly control (all inputs at 6.62x).

    Answers "how close does each simulated assembly get to the real ONT
    assembly", which the 10x main table cannot answer because the real ONT
    control only carries 6.62x.
    """
    dirs = table_dirs(project)
    labels = {
        "real_ont": "Real ONT",
        "errhmm": "GCerrHMM",
        "nanosim": "NanoSim",
        "pbsim": "PBSim3-sample",
        "badread": "badread",
    }
    rows = []
    reference = None
    for tool, label in labels.items():
        path = next(
            (
                hit
                for directory in dirs
                if (
                    hit := find_table(
                        [directory],
                        f"assembly_flye_{tool}_Hsapiens_chr21_6p6x.json",
                    )
                )
            ),
            None,
        )
        if path is None:
            continue
        payload = read_json(path)
        row = {
            "species": "Hsapiens_chr21",
            "tool_key": tool,
            "tool_label": label,
            "source": path.name,
            "num_contigs": payload.get("num_contigs", ""),
            "total_length": payload.get("total_length", ""),
            "n50": payload.get("n50", ""),
            "reference_identity": payload.get("reference_identity", ""),
        }
        if label == "Real ONT":
            reference = row
        rows.append(row)

    if reference:
        for row in rows:
            for metric in ("reference_identity", "n50", "num_contigs"):
                value = row.get(metric)
                base = reference.get(metric)
                if isinstance(value, (int, float)) and isinstance(base, (int, float)):
                    row[f"delta_{metric}_vs_real"] = round(value - base, 4)
                else:
                    row[f"delta_{metric}_vs_real"] = ""

    fields = [
        "species", "tool_key", "tool_label", "source",
        "num_contigs", "total_length", "n50", "reference_identity",
        "delta_reference_identity_vs_real",
        "delta_n50_vs_real",
        "delta_num_contigs_vs_real",
    ]
    write_csv(project / "results/framework/stats/assembly_coverage_matched.csv", fields, rows)
    return rows


def build_sv_table(project):
    """SV spike-in results, one row per species x tool.

    Multiple spike-in designs (seeds) are averaged; the per-design values are
    kept in the row so the paper can report dispersion instead of a single
    point estimate.
    """
    dirs = table_dirs(project)
    rows = []
    for species in EXTERNAL_SPECIES:
        for tool, label in EXTERNAL_TOOLS.items():
            found = []
            for directory in dirs:
                found = sorted(
                    directory.glob(f"variant_sv_{tool}_{species}_spikein*.json")
                )
                if found:
                    break
            if not found:
                continue
            payloads = [(p.name, read_json(p)) for p in found]
            overalls = [payload.get("overall", {}) for _, payload in payloads]
            row = {
                "species": species,
                "tool_key": tool,
                "tool_label": label,
                "source": ";".join(name for name, _ in payloads),
                "n_designs": len(payloads),
                "sv_caller": payloads[0][1].get("sv_caller", ""),
                "min_svlen": payloads[0][1].get("min_svlen", ""),
                # Designs can carry different event counts (the 39-event pilot
                # alongside the 104-event replicates), so the pooled truth is
                # the sum, and the per-design counts stay visible.
                "truth_total": sum(
                    int(payload.get("truth_total") or 0) for _, payload in payloads
                ),
                "truth_per_design": ",".join(
                    str(payload.get("truth_total", "")) for _, payload in payloads
                ),
            }
            for metric in ("precision", "recall", "f1"):
                values = [o.get(metric) for o in overalls if o.get(metric) is not None]
                mean, sd = mean_sd(values)
                key = "f1_score" if metric == "f1" else metric
                row[key] = round(mean, 4) if values else ""
                row[f"{key}_sd"] = round(sd, 4) if len(values) > 1 else 0.0
            row["called_total"] = (
                round(mean_sd([p.get("called_total") for _, p in payloads])[0], 1)
                if payloads
                else ""
            )
            tp = sum(o.get("tp", 0) for o in overalls)
            fp = sum(o.get("fp", 0) for o in overalls)
            fn = sum(o.get("fn", 0) for o in overalls)
            row["pooled_tp"], row["pooled_fp"], row["pooled_fn"] = tp, fp, fn
            row["pooled_recall_ci95"] = wilson_ci(tp, tp + fn)
            row["pooled_precision_ci95"] = wilson_ci(tp, tp + fp)

            type_names = sorted(
                {
                    svtype
                    for _, payload in payloads
                    for svtype in payload.get("by_type", {})
                }
            )
            for svtype in type_names:
                for metric in ("precision", "recall", "f1"):
                    values = [
                        payload["by_type"][svtype].get(metric)
                        for _, payload in payloads
                        if svtype in payload.get("by_type", {})
                    ]
                    mean, sd = mean_sd(values)
                    row[f"{svtype.lower()}_{metric}"] = (
                        round(mean, 4) if values else ""
                    )
                    row[f"{svtype.lower()}_{metric}_sd"] = (
                        round(sd, 4) if len(values) > 1 else 0.0
                    )
            label_bins = sorted(
                {
                    label_bin
                    for _, payload in payloads
                    for label_bin in payload.get("by_size", {})
                }
            )
            for label_bin in label_bins:
                values = [
                    payload["by_size"][label_bin].get("f1")
                    for _, payload in payloads
                    if label_bin in payload.get("by_size", {})
                ]
                mean, sd = mean_sd(values)
                key = "size_" + label_bin.replace("-", "_").replace("+", "plus")
                row[f"{key}_f1"] = round(mean, 4) if values else ""
                row[f"{key}_f1_sd"] = round(sd, 4) if len(values) > 1 else 0.0
            rows.append(row)
    extra = sorted(
        {
            key
            for row in rows
            for key in row
            if key not in {
                "species", "tool_key", "tool_label", "source", "sv_caller",
                "min_svlen", "truth_total", "truth_per_design",
                "called_total", "n_designs",
                "precision", "recall", "f1_score",
                "precision_sd", "recall_sd", "f1_score_sd",
                "pooled_tp", "pooled_fp", "pooled_fn",
                "pooled_recall_ci95", "pooled_precision_ci95",
            }
        }
    )
    fields = [
        "species", "tool_key", "tool_label", "source", "sv_caller",
        "min_svlen", "truth_total", "truth_per_design", "called_total", "n_designs",
        "precision", "precision_sd", "recall", "recall_sd",
        "f1_score", "f1_score_sd",
        "pooled_tp", "pooled_fp", "pooled_fn",
        "pooled_recall_ci95", "pooled_precision_ci95",
    ] + extra
    write_csv(project / "results/framework/stats/sv_calling.csv", fields, rows)
    return rows


def build_hifi_table(project):
    grouped = {}
    for path in sorted(project.glob("results/gc_improvement_hifi/**/gc_improvement_summary.json")):
        payload = read_json(path)
        metadata = payload.get("metadata", {})
        if float(metadata.get("coverage", 0)) != 10.0:
            continue
        species = metadata.get("species", path.parents[1].name)
        for route in payload.get("routes", []):
            key = (species, route.get("route", ""))
            bucket = grouped.setdefault(key, {
                "species": species,
                "platform": metadata.get("platform", "hifi"),
                "coverage": metadata.get("coverage", 10.0),
                "route": route.get("route", ""),
                "composite": [],
                "read_length": [],
                "qv": [],
                "gc": [],
                "kmer": [],
                "source": path.name,
            })
            for metric in ("composite", "read_length", "qv", "gc", "kmer"):
                bucket[metric].append(route.get(metric))
    rows = []
    for bucket in grouped.values():
        row = {
            "species": bucket["species"],
            "platform": bucket["platform"],
            "coverage": bucket["coverage"],
            "route": bucket["route"],
            "source": bucket["source"],
        }
        for metric in ("composite", "read_length", "qv", "gc", "kmer"):
            mean, sd = mean_sd(bucket[metric])
            row[f"{metric}_mean"] = mean
            row[f"{metric}_sd"] = sd
        row["n"] = len(bucket["composite"])
        rows.append(row)
    fields = [
        "species", "platform", "coverage", "route", "n",
        "composite_mean", "composite_sd",
        "read_length_mean", "read_length_sd",
        "qv_mean", "qv_sd",
        "gc_mean", "gc_sd",
        "kmer_mean", "kmer_sd",
        "source",
    ]
    write_csv(project / "results/framework/stats/hifi_level1_10x.csv", fields, rows)
    return rows


def build_status_table(project, reads, mapping, variant, assembly, sv):
    keys = set()
    for rows in (reads, mapping, variant, sv, assembly):
        keys.update((row["species"], row["tool_label"]) for row in rows)
    status_rows = []
    for species, tool_label in sorted(keys):
        status_rows.append({
            "species": species,
            "tool_label": tool_label,
            "reads": "measured" if any(r["species"] == species and r["tool_label"] == tool_label for r in reads) else "missing",
            "mapping": "measured" if any(r["species"] == species and r["tool_label"] == tool_label for r in mapping) else "missing",
            "variant": "measured" if any(r["species"] == species and r["tool_label"] == tool_label for r in variant) else "missing",
            "sv": "measured" if any(r["species"] == species and r["tool_label"] == tool_label for r in sv) else "missing",
            "assembly": "measured" if any(r["species"] == species and r["tool_label"] == tool_label for r in assembly) else "missing",
        })
    fields = ["species", "tool_label", "reads", "mapping", "variant", "sv", "assembly"]
    write_csv(project / "results/framework/stats/four_layer_status.csv", fields, status_rows)
    return status_rows


def write_summary(project, reads, mapping, variant, sv, assembly, status, hifi):
    path = project / "results/framework/stats/FRAMEWORK_SUMMARY.md"
    lines = [
        "# Teacher framework summary",
        "",
        "| Layer | Rows | External species panel |",
        "|---|---:|---|",
        f"| Reads / Level-1 | {len(reads)} | Ecoli, Scerevisiae, Hsapiens_chr21 |",
        f"| Mapping | {len(mapping)} | external rows appear as they complete |",
        f"| Variant calling | {len(variant)} | external rows appear as they complete |",
        f"| SV spike-in calling | {len(sv)} | Ecoli, Scerevisiae (donor with injected SVs) |",
        f"| Assembly | {len(assembly)} | external rows appear as they complete |",
        f"| HiFi cross-platform | {len(hifi)} | GCerrHMM route comparison |",
        "",
        "## Layer status",
        "",
        "| Species | Tool | Reads | Mapping | Variant | SV | Assembly |",
        "|---|---|---|---|---|---|---|",
    ]
    for row in status:
        lines.append(
            f"| {row['species']} | {row['tool_label']} | {row['reads']} | "
            f"{row['mapping']} | {row['variant']} | {row['sv']} | {row['assembly']} |"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"FRAMEWORK_SUMMARY={path}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()
    project = Path(args.project_dir).resolve()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    reads = build_reads_table(project)
    mapping = build_mapping_table(project)
    variant = build_variant_table(project)
    sv = build_sv_table(project)
    assembly = build_assembly_table(project)
    coverage_matched = build_coverage_matched_table(project)
    hifi = build_hifi_table(project)
    status = build_status_table(project, reads, mapping, variant, assembly, sv)
    write_summary(project, reads, mapping, variant, sv, assembly, status, hifi)
    print(f"COVERAGE_MATCHED_ROWS={len(coverage_matched)}")
    print(f"FRAMEWORK_TABLES={output_dir}")


if __name__ == "__main__":
    main()
