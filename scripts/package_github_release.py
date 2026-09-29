#!/usr/bin/env python3
"""Assemble a clean, portable, GitHub-ready GCerrHMM release directory."""

from __future__ import annotations

import argparse
import re
import shutil
import zipfile
from pathlib import Path

from common_io import sha256_file


ROOT_FILES = [
    ".gitattributes",
    ".gitignore",
    ".zenodo.json",
    "CITATION.cff",
    "LICENSE",
    "README.md",
    "REPRODUCIBILITY.md",
    "UPLOAD_GUIDE.md",
    "environment-tools.yml",
    "environment.yml",
    "reproduce.sh",
    "requirements.txt",
]

OPTIONAL_ROOT_FILES = [
    "CHANGELOG.md",
    "CODE_OF_CONDUCT.md",
    "CONTRIBUTING.md",
    "Makefile",
    "pyproject.toml",
]

CORE_SCRIPTS = [
    "ablation_claim_decision.py",
    "apply_author_metadata.py",
    "audit_code_similarity.py",
    "audit_data_integrity.py",
    "audit_hardcoded_paths.py",
    "audit_logs.py",
    "audit_panel_provenance.py",
    "audit_w3_submission.py",
    "bai_ref_span.py",
    "build_additional_files.py",
    "build_background_figures.py",
    "build_chinese_abstract.py",
    "build_chinese_full_manuscript.py",
    "build_delta_to_real_figure.py",
    "build_figure_contact_sheet.py",
    "build_gcerrhmm_main_figures.py",
    "build_r1_level1_overall.py",
    "build_teacher_framework_figures.py",
    "build_v2_layer_summary.py",
    "build_w3_review_copy.py",
    "build_w3_submission.py",
    "check_gc_model.py",
    "common_io.py",
    "coverage_power_analysis.py",
    "delta_to_real_decision.py",
    "diagnose_kmer_overlap.py",
    "docx_style.py",
    "download_fastq_multipart.py",
    "download_verified_sources.py",
    "error_spectrum.py",
    "evaluate_assembly_fasta.py",
    "export_framework_tables.py",
    "export_gc_bins_cross_table.py",
    "figure_io.py",
    "gc_error_fidelity.py",
    "innovation_postanalysis.py",
    "package_github_release.py",
    "plot_gc_error_curve.py",
    "pooled_gc_curve.py",
    "project_paths.py",
    "recompute_kmer_metric.py",
    "reproduce_gc_improvement.py",
    "retrain_species.py",
    "select_mapq_threshold.py",
    "static_audit.py",
    "subsample_fastq_to_bases.py",
    "summarize_hifi_seeds.py",
    "verify_freeze_manifest.py",
]

DOCS = [
    "BMC_Technical_Design_v2.md",
    "claim_wording_options.md",
    "composite_evaluation_protocol.md",
    "cross_tool_reproducibility_notes.md",
    "delta_to_real_framing.md",
    "data_freeze_manifest.md",
    "errhmm_innovation_figures_20260918.zip",
    "errhmm_paper_visuals_20260918.zip",
    "figure_captions.md",
    "final_cross_tool_comparison.md",
    "framework_data_gaps.md",
    "gc_aware_algorithm.md",
    "GITHUB_RELEASE_CHECKLIST.md",
    "gcerrhmm_main_figure_captions.md",
    "gcerrhmm_main_figure_status.md",
    "hifi_region_strategy.md",
    "kmer_metric_audit.md",
    "panel_provenance.md",
    "REVIEWER_GUIDE.md",
    "REPRODUCIBILITY_AUDIT_20260927.md",
    "sv_evaluation_protocol.md",
    "teacher_framework_design.md",
]

REPRODUCIBILITY_FILES = [
    "checksums.sha256",
    "gc_improvement_summary.csv",
    "gc_improvement_summary.json",
    "gc_improvement_summary.md",
    "r1_level1_overall.csv",
    "r1_level1_overall.md",
    "r1_level1_overall.sha256",
    "r1_provenance.md",
]

MANUSCRIPT_FILES = [
    "AUTHOR_INPUT_REQUIRED.md",
    "GCerrHMM_BMC_Research_article_W3.3_review_with_figures_20260929.pdf",
    "GCerrHMM_BMC_Research_article_W3.3_submission_20260929.pdf",
    "GCerrHMM_Chinese_figure_abstract_W3.3_20260929.pdf",
    "GCerrHMM_Chinese_full_manuscript_W3.3_20260929.pdf",
    "W3.2_review_response.md",
    "W3.3_figures_contact_sheet.png",
    "W3_reference_audit.md",
    "W3_style_notes.md",
    "author_metadata.template.json",
]

TEXT_SUFFIXES = {
    ".cff",
    ".cfg",
    ".csv",
    ".in",
    ".json",
    ".md",
    ".py",
    ".sh",
    ".toml",
    ".txt",
    ".yaml",
    ".yml",
}

PORTABLE_REPLACEMENTS = (
    (re.compile(r"/mnt/c/Users/[^/]+/Desktop/[^/]+/?"), ""),
    (re.compile(r"[A-Za-z]:[\\/]Users[\\/][^\\/]+[\\/]Desktop[\\/][^\\/]+[\\/]?"), ""),
    (re.compile(r"/home/[^/\s]+/errhmm_project/?"), ""),
    (re.compile(r"/home/[^/\s]+/bmc_data/?"), "<external-data>/"),
    (re.compile(r"/home/[^/\s]+/miniconda3"), "$CONDA_PREFIX"),
    (re.compile(r"/home/[^/\s]+/apt-bcftools/root/usr/bin"), "<bcftools-prefix>"),
    (re.compile(r"\b[\w.-]+@10\.70\.5\.64\b"), "<remote-host>"),
    (re.compile(r"\b10\.70\.5\.64\b"), "<remote-host>"),
    (re.compile(r"\bbmc_project\b"), "project"),
)

PRIVATE_USER_MARKERS = (
    b"27" + b"947",
    b"bob" + b"by",
    b"z" + b"by",
)

PRIVATE_PATH_MARKERS = (
    b"27" + b"947",
    (b"10." + b"70.5.64"),
    b"/home/" + b"bob" + b"by",
    b"C:\\Users\\" + b"27" + b"947",
    b"/mnt/c/Users/" + b"27" + b"947",
    b"Desktop/bmc_" + b"project",
)


def copy_portable(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if source.suffix.lower() not in TEXT_SUFFIXES:
        shutil.copy2(source, destination)
        return
    try:
        text = source.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        # Some frozen CSV/Markdown artifacts are legacy GBK. Preserve their
        # bytes (and therefore their recorded checksums) and audit raw bytes.
        shutil.copy2(source, destination)
        return
    for pattern, replacement in PORTABLE_REPLACEMENTS:
        text = pattern.sub(replacement, text)
    destination.write_bytes(text.encode("utf-8"))


def ensure_safe_output(source: Path, output: Path) -> None:
    if output == source or source not in output.parents:
        raise SystemExit(
            f"refusing to replace output outside source tree: {output}"
        )


def scan_private_markers(root: Path) -> list[str]:
    findings = []
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        payload = path.read_bytes()
        markers = (
            PRIVATE_USER_MARKERS
            if path.suffix.lower() in TEXT_SUFFIXES
            else PRIVATE_PATH_MARKERS
        )
        for marker in markers:
            if marker in payload:
                findings.append(
                    f"{path.relative_to(root).as_posix()}: {marker.decode()}"
                )
                break
    return findings


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", default=".")
    parser.add_argument("--output", required=True)
    parser.add_argument("--zip", dest="zip_path")
    args = parser.parse_args()

    source = Path(args.source).resolve()
    output = Path(args.output).resolve()
    ensure_safe_output(source, output)

    staging = output.with_name(output.name + ".building")
    if staging.exists():
        shutil.rmtree(staging)
    staging.mkdir(parents=True)

    for filename in ROOT_FILES:
        copy_portable(source / filename, staging / filename)
    for filename in OPTIONAL_ROOT_FILES:
        path = source / filename
        if path.exists():
            copy_portable(path, staging / filename)

    workflow = source / ".github" / "workflows" / "smoke.yml"
    copy_portable(workflow, staging / ".github" / "workflows" / "smoke.yml")

    for path in sorted((source / "src").glob("*.py")):
        copy_portable(path, staging / "src" / path.name)
    for path in sorted((source / "tests").glob("*.py")):
        copy_portable(path, staging / "tests" / path.name)
    for path in sorted((source / "config").glob("*.json")):
        copy_portable(path, staging / "config" / path.name)
    for path in sorted((source / "config" / "read_profiles").glob("*.json")):
        copy_portable(
            path,
            staging / "config" / "read_profiles" / path.name,
        )
    for filename in CORE_SCRIPTS:
        copy_portable(
            source / "scripts" / filename,
            staging / "scripts" / filename,
        )

    for filename in DOCS:
        copy_portable(source / "docs" / filename, staging / "docs" / filename)
    for filename in REPRODUCIBILITY_FILES:
        copy_portable(
            source / "docs" / "reproducibility" / filename,
            staging / "docs" / "reproducibility" / filename,
        )
    for filename in MANUSCRIPT_FILES:
        copy_portable(
            source / "docs" / "manuscript" / filename,
            staging / "docs" / "manuscript" / filename,
        )

    for directory in (
        "docs/background_figures",
        "docs/gcerrhmm_main_figures_refined_20260927_v2",
    ):
        for path in sorted((source / directory).glob("*")):
            if path.is_file():
                copy_portable(path, staging / directory / path.name)

    for path in sorted(staging.rglob("__pycache__"), reverse=True):
        if path.is_dir():
            shutil.rmtree(path)

    findings = scan_private_markers(staging)
    if findings:
        shutil.rmtree(staging)
        raise SystemExit(
            "private path audit failed:\n- " + "\n- ".join(findings)
        )

    manifest_lines = []
    for path in sorted(staging.rglob("*")):
        if path.is_file() and path.name != "RELEASE_MANIFEST.sha256":
            relative = path.relative_to(staging)
            manifest_lines.append(f"{sha256_file(path)}  {relative.as_posix()}")
    (staging / "RELEASE_MANIFEST.sha256").write_text(
        "\n".join(manifest_lines) + "\n",
        encoding="utf-8",
    )

    if output.exists():
        shutil.rmtree(output)
    staging.replace(output)

    if args.zip_path:
        zip_path = Path(args.zip_path).resolve()
        zip_path.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as archive:
            for path in sorted(output.rglob("*")):
                if path.is_file():
                    info = zipfile.ZipInfo(
                        path.relative_to(output).as_posix(),
                        date_time=(2026, 9, 29, 0, 0, 0),
                    )
                    info.compress_type = zipfile.ZIP_DEFLATED
                    info.external_attr = 0o100644 << 16
                    archive.writestr(info, path.read_bytes())
        print(f"ZIP={zip_path}")
    print(f"RELEASE={output}")


if __name__ == "__main__":
    main()
