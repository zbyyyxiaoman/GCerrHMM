#!/usr/bin/env python3
"""Build the frozen six-species x three-route Level-1 overall table."""

import argparse
import csv
import glob
import hashlib
import json
import statistics
from pathlib import Path


SPECIES = [
    ("Ecoli", "E. coli"),
    ("Scerevisiae", "S. cerevisiae"),
    ("Athaliana", "A. thaliana"),
    ("Dmelanogaster", "D. melanogaster"),
    ("Mmusculus_chr19", "M. musculus chr19"),
    ("Hsapiens_chr21", "H. sapiens chr21"),
]
ROUTES = [
    ("route_A_sample", "A: Sample"),
    ("route_B_qshmm", "B: qsHMM (1-bin)"),
    ("route_C_errhmm", "C: errHMM (GC-aware)"),
]
WEIGHTS = {
    "read_length": 0.15,
    "qv": 0.20,
    "gc": 0.25,
    "kmer": 0.15,
}


def last_line_number(lines, needle):
    matches = [index + 1 for index, line in enumerate(lines) if needle in line]
    return matches[-1] if matches else None


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", required=True)
    args = parser.parse_args()
    project_dir = Path(args.project_dir)
    tables_dir = project_dir / "results" / "tables"
    stats_dir = project_dir / "results" / "stats"
    stats_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    provenance = []
    all_inputs = []
    max_delta = 0.0
    for species, species_label in SPECIES:
        for route, route_label in ROUTES:
            paths = sorted(
                (tables_dir / name)
                for name in glob.glob(
                    str(tables_dir / f"level1_{species}_{route}_r*.json")
                )
            )
            if not paths:
                raise SystemExit(f"Missing frozen result for {species}/{route}")
            values = []
            kmer_values = []
            source_files = []
            for path in paths:
                raw_lines = path.read_text().splitlines()
                data = json.loads(path.read_text())
                composite = data["composite_score"]
                sub_scores = composite["sub_scores"]
                file_score = float(composite["composite_score"])
                locked_score = sum(
                    WEIGHTS[metric] * float(sub_scores[metric])
                    for metric in WEIGHTS
                ) / sum(WEIGHTS.values())
                delta = abs(file_score - locked_score)
                max_delta = max(max_delta, delta)
                values.append(locked_score)
                if "kmer" in sub_scores:
                    kmer_values.append(float(sub_scores["kmer"]))
                source_files.append(str(path.relative_to(project_dir)))
                all_inputs.append(path)
                provenance.append({
                    "species": species_label,
                    "route": route_label,
                    "replicate": path.stem.rsplit("_r", 1)[-1],
                    "source_file": str(path.relative_to(project_dir)),
                    "composite_line": last_line_number(
                        raw_lines, '"composite_score": {'
                    ),
                    "sub_scores_line": last_line_number(
                        raw_lines, '"sub_scores": {'
                    ),
                    "file_score": file_score,
                    "locked_score": locked_score,
                    "delta": delta,
                })
            mean = statistics.mean(values)
            sd = statistics.stdev(values) if len(values) > 1 else 0.0
            kmer_spread = (max(kmer_values) - min(kmer_values)
                           if len(kmer_values) > 1 else 0.0)
            kmer_sd = (statistics.stdev(kmer_values)
                       if len(kmer_values) > 1 else 0.0)
            rows.append({
                "species": species_label,
                "route": route_label,
                "n": len(values),
                "seeds_or_replicates": "legacy r1/r2; seed not persisted "
                "(generator default seed=42)",
                "mean": mean,
                "sd": sd,
                "raw_values": ";".join(f"{value:.6f}" for value in values),
                "kmer_mean": (statistics.mean(kmer_values)
                              if kmer_values else None),
                "kmer_spread": kmer_spread,
                "kmer_sd": kmer_sd,
                "source_files": ";".join(source_files),
                "formula": "normalized locked weights "
                "0.15 read_length + 0.20 qv + 0.25 gc + 0.15 kmer; "
                "error-rate disabled",
            })

    csv_path = stats_dir / "r1_level1_overall.csv"
    md_path = stats_dir / "r1_level1_overall.md"
    provenance_path = stats_dir / "r1_provenance.md"

    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "species", "route", "n", "seeds_or_replicates",
                "mean", "sd", "raw_values", "kmer_mean", "kmer_spread",
                "kmer_sd", "source_files", "formula",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    with md_path.open("w") as handle:
        handle.write("# R1 Level-1 overall\n\n")
        handle.write(
            "Frozen aggregation only. Composite uses the locked weights "
            "with error-rate disabled and renormalization.\n\n"
        )
        handle.write(
            "| Species | Route | n | Seeds/replicates | Mean ± SD | Raw values "
            "| k-mer mean | k-mer replicate spread |\n"
        )
        handle.write("|---|---|---:|---|---:|---|---:|---:|\n")
        for row in rows:
            kmer_mean = ("-" if row["kmer_mean"] is None
                         else f"{row['kmer_mean']:.2f}")
            handle.write(
                f"| {row['species']} | {row['route']} | {row['n']} | "
                f"{row['seeds_or_replicates']} | "
                f"{row['mean']:.4f} ± {row['sd']:.4f} | "
                f"{row['raw_values']} | {kmer_mean} | "
                f"{row['kmer_spread']:.2f} |\n"
            )
        bc_spread = max((r["kmer_spread"] for r in rows
                         if r["route"].startswith(("B:", "C:"))), default=0.0)
        any_spread = max((r["kmer_spread"] for r in rows), default=0.0)
        handle.write(
            "\n**Read this table with the k-mer column open.** The composite "
            "contains no term that measures the error rate conditional on local "
            "GC, and its k-mer term is noisy: the replicate-to-replicate spread "
            "of that sub-score reaches "
            f"{bc_spread:.1f} points on a 0-100 scale within a single model "
            f"route (B/C) and {any_spread:.1f} points on the sampling baseline. "
            "Comparing the GC-aware route with the "
            "1-bin control on this table is therefore reported as **not "
            "discriminable**, not as a win for either route; see "
            "`docs/gc_claim_evidence_chain.md` for the pre-registered "
            "instruments that do carry that claim.\n"
        )

    with provenance_path.open("w") as handle:
        handle.write("# R1 provenance\n\n")
        handle.write(
            "Every score is recomputed from the frozen Level-1 sub-score "
            "fields in the listed JSON file. The normalized locked weights "
            "are `0.15 read_length + 0.20 qv + 0.25 gc + 0.15 kmer`.\n\n"
        )
        handle.write(
            "| Species | Route | Replicate | Frozen file | Composite line | "
            "Sub-scores line | File composite | Locked recompute | Delta |\n"
        )
        handle.write("|---|---|---:|---|---:|---:|---:|---:|---:|\n")
        for row in provenance:
            handle.write(
                f"| {row['species']} | {row['route']} | {row['replicate']} | "
                f"`{row['source_file']}` | {row['composite_line']} | "
                f"{row['sub_scores_line']} | {row['file_score']:.6f} | "
                f"{row['locked_score']:.6f} | {row['delta']:.6f} |\n"
            )
        handle.write(
            f"\nMaximum file-vs-recompute delta: `{max_delta:.8f}`. "
            "The same sub-score fields are used by fig3; the k-mer "
            "C-minus-B values used by fig4 are recomputed from the same "
            "frozen `kmer` sub-scores.\n"
        )
        handle.write(
            "\n**Paired k-mer bootstrap.** `results/stats/kmer_pair_summary.csv` "
            "re-estimates the k-mer correlation with 1000 reads per side and a "
            "paired bootstrap over reads (200 resamples, same read indices for "
            "both routes). Eleven of twelve route/replicate comparisons have a "
            "difference interval containing zero; the single interval that "
            "excludes zero favours the 1-bin control and does not reproduce in "
            "the second replicate of the same genome. The six-species panel is "
            "therefore reported as not discriminable on this statistic.\n"
        )

    checksum_path = stats_dir / "r1_level1_overall.sha256"
    checksum_paths = [csv_path, md_path, provenance_path, *all_inputs]
    checksum_lines = [
        f"{sha256(path)}  {path.relative_to(project_dir)}"
        for path in checksum_paths
    ]
    checksum_path.write_text("\n".join(checksum_lines) + "\n")

    # The frozen manifest lives in the code checkout; when the analysis runs
    # against a results-only project tree, fall back to the repository copy
    # instead of crashing at the last step.
    candidates = [
        project_dir / "docs" / "data_freeze_manifest.md",
        Path(__file__).resolve().parents[1] / "docs" / "data_freeze_manifest.md",
    ]
    manifest_path = next((p for p in candidates if p.exists()), None)
    if manifest_path is None:
        print("WARNING: data_freeze_manifest.md not found; skipping append")
        return
    marker = "<!-- R1_LEVEL1_OVERALL_SHA256 -->"
    manifest_text = manifest_path.read_text()
    if marker not in manifest_text:
        manifest_path.write_text(
            manifest_text.rstrip()
            + "\n\n"
            + marker
            + "\n## R1 Level-1 overall SHA-256\n\n"
            + "\n".join(f"- `{line}`" for line in checksum_lines)
            + "\n"
        )

    print(f"R1_CSV={csv_path}")
    print(f"R1_MD={md_path}")
    print(f"R1_PROVENANCE={provenance_path}")
    print(f"R1_SHA256={checksum_path}")
    print(f"R1_MAX_DELTA={max_delta:.8f}")


if __name__ == "__main__":
    main()
