#!/usr/bin/env python3
"""Build the Additional file 1-7 package for the BMC submission."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.shared import Inches, Pt
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill
from openpyxl.utils import get_column_letter


def read_csv(path: Path) -> list[list[str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.reader(handle))


def find_input(project: Path, package: Path, relative: Path) -> Path:
    candidates = (
        package / relative,
        project / relative,
    )
    for candidate in candidates:
        if candidate.exists() and candidate.stat().st_size > 0:
            return candidate
    raise SystemExit(f"missing additional-file input: {relative}")


def write_workbook(
    output: Path,
    sheets: list[tuple[str, Path]],
) -> None:
    workbook = Workbook()
    workbook.remove(workbook.active)
    for sheet_name, source in sheets:
        rows = read_csv(source)
        sheet = workbook.create_sheet(sheet_name)
        for row in rows:
            sheet.append(row)
        for cell in sheet[1]:
            cell.font = Font(bold=True, color="FFFFFF")
            cell.fill = PatternFill("solid", fgColor="1F4E79")
        for column_index in range(1, len(rows[0]) + 1):
            letter = get_column_letter(column_index)
            width = max(
                len(str(row[column_index - 1]))
                for row in rows
            ) + 2
            sheet.column_dimensions[letter].width = min(max(width, 10), 42)
        sheet.freeze_panes = "A2"
    output.parent.mkdir(parents=True, exist_ok=True)
    workbook.save(output)


def add_paragraph(
    document: Document,
    text: str,
    style: str | None = None,
    bold: bool = False,
) -> None:
    paragraph = document.add_paragraph(style=style)
    run = paragraph.add_run(text)
    run.bold = bold
    return paragraph


def set_note_style(document: Document) -> None:
    normal = document.styles["Normal"]
    normal.font.name = "Arial"
    normal.font.size = Pt(10.5)
    normal.paragraph_format.space_after = Pt(6)
    for name, size in (
        ("Title", 16),
        ("Heading 1", 13),
        ("Heading 2", 11),
    ):
        try:
            style = document.styles[name]
        except KeyError:
            continue
        style.font.name = "Arial"
        style.font.size = Pt(size)
        style.font.bold = True


def write_note(
    output: Path,
    title: str,
    sections: list[tuple[str, list[str]]],
    table: list[list[str]] | None = None,
) -> None:
    document = Document()
    set_note_style(document)
    heading = document.add_paragraph(style="Title")
    heading.add_run(title)
    for section_title, paragraphs in sections:
        document.add_heading(section_title, level=1)
        for text in paragraphs:
            add_paragraph(document, text)
    if table:
        document.add_heading("Data used", level=1)
        table_object = document.add_table(rows=1, cols=len(table[0]))
        table_object.style = "Table Grid"
        for index, value in enumerate(table[0]):
            table_object.rows[0].cells[index].text = str(value)
        for row in table[1:]:
            cells = table_object.add_row().cells
            for index, value in enumerate(row):
                cells[index].text = str(value)
    output.parent.mkdir(parents=True, exist_ok=True)
    document.save(output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", required=True)
    parser.add_argument("--package-dir", default=None)
    parser.add_argument("--output-dir", default=None)
    args = parser.parse_args()

    project = Path(args.project_dir)
    package = Path(args.package_dir) if args.package_dir else project / "docs"
    output = (
        Path(args.output_dir)
        if args.output_dir
        else project / "docs" / "additional_files"
    )
    output.mkdir(parents=True, exist_ok=True)

    write_workbook(
        output / "Additional_file_1_species_panel.xlsx",
        [(
            "Species panel",
            find_input(
                project,
                package,
                Path("docs/paper_figures/table1_species_panel.csv"),
            ),
        )],
    )
    write_workbook(
        output / "Additional_file_2_level1_sub_scores.xlsx",
        [
            (
                "ONT_frozen",
                find_input(
                    project,
                    package,
                    Path("results/framework/stats/reads_level1.csv"),
                ),
            ),
            (
                "HiFi_10x",
                find_input(
                    project,
                    package,
                    Path("results/framework/stats/hifi_level1_10x.csv"),
                ),
            ),
            (
                "HiFi_seed_deltas",
                find_input(
                    project,
                    package,
                    Path("results/framework/stats/hifi_crossplatform_seeds.csv"),
                ),
            ),
        ],
    )
    write_workbook(
        output / "Additional_file_3_gc_bin_ablation.xlsx",
        [
            (
                "Decision",
                find_input(
                    project,
                    package,
                    Path("results/stats/gc_bins_claim_decision.csv"),
                ),
            ),
            (
                "Seed_details",
                find_input(
                    project,
                    package,
                    Path(
                        "results/stats/"
                        "profile_matched_gc_bins_seed_details_20260916_010047.csv"
                    ),
                ),
            ),
        ],
    )
    write_workbook(
        output / "Additional_file_4_delta_to_real.xlsx",
        [
            (
                "Decision",
                find_input(
                    project,
                    package,
                    Path(
                        "results/framework/stats/delta_to_real_decision.csv"
                    ),
                ),
            ),
            (
                "Panel_30x",
                find_input(
                    project,
                    package,
                    Path(
                        "results/framework/stats/delta_to_real_panel_30x.csv"
                    ),
                ),
            ),
        ],
    )

    write_note(
        output / "Additional_file_5_instrument_audit.docx",
        "Additional file 5. Instrument audit: four measurement problems, "
        "symptoms and fixes",
        [
            (
                "Purpose",
                [
                    "This note records four measurement problems that produced "
                    "plausible but wrong numbers during development. Each was "
                    "fixed in code and is now protected by an executable audit "
                    "or fail-loud gate.",
                ],
            ),
            (
                "Problem 1: MAPQ threshold",
                [
                    "A default minimum MAPQ of 20 removed reads carrying a large "
                    "share of the chr21 error signal. The 30x chr21 model was "
                    "retrained with minimum MAPQ 0, and the threshold is now "
                    "reported explicitly.",
                ],
            ),
            (
                "Problem 2: insertion-length cap",
                [
                    "A single insertion state made every insertion one base long. "
                    "Run-length states I1-I3/I4+ and D1-D3/D4+ remove this "
                    "structural cap and are used for the reported model.",
                ],
            ),
            (
                "Problem 3: silent profile fallback",
                [
                    "A missing read-length profile silently fell back to a "
                    "platform default. The entry gate now refuses to start "
                    "without a read-length and quality profile, and the exit "
                    "gate verifies produced read count, total bases and mean "
                    "read length within 10% of the request.",
                ],
            ),
            (
                "Problem 4: k-mer estimator spread",
                [
                    "At 250 reads per side the k-mer sub-score had a replicate "
                    "spread larger than the between-route differences it was used "
                    "to rank. The reported estimator uses 1000 reads per side "
                    "and a paired bootstrap over reads; the k-mer term is treated "
                    "as descriptive rather than as the primary GC claim.",
                ],
            ),
            (
                "Verification",
                [
                    "The verification entry points are `reproduce.sh --tests`, "
                    "`reproduce.sh --static-check` and `reproduce.sh --smoke`; "
                    "the reviewer checklist is docs/REVIEWER_GUIDE.md.",
                ],
            ),
        ],
    )

    write_note(
        output / "Additional_file_6_uncertainty_layers.docx",
        "Additional file 6. Two layers of uncertainty in the GC-fidelity "
        "measurement",
        [
            (
                "Read sampling",
                [
                    "Disjoint subsamples of a pooled 30x simulation gave a "
                    "correlation SD of 0.176 at 5x, 0.194 at 10x and 0.009 at "
                    "15x. This component answers how deeply the reads must be "
                    "sampled and is the only component interpreted as a power "
                    "statement.",
                ],
            ),
            (
                "Between-window heterogeneity",
                [
                    "The between-window bootstrap width was 1.01, 0.96 and 0.96 "
                    "at 10x, 20x and 30x, respectively. It did not shrink with "
                    "depth and is not a power statement. Conflating it with read "
                    "sampling makes the same data appear both underpowered and "
                    "precise.",
                ],
            ),
            (
                "Interpretation",
                [
                    "At 15x the read-sampling component is small enough to "
                    "exclude an effect above approximately 0.05 in correlation "
                    "units. The reported GC effect is therefore an effect-size "
                    "statement at 15x and above, not a claim that the experiment "
                    "lacks precision at every depth.",
                ],
            ),
        ],
    )

    write_note(
        output / "Additional_file_7_hifi_crossplatform_note.docx",
        "Additional file 7. HiFi cross-platform chain and the "
        "encoding-limited QV caveat",
        [
            (
                "Cross-platform chain",
                [
                    "The framework supports an ONT primary chain and a PacBio "
                    "HiFi cross-platform chain. HiFi assembly uses hifiasm, "
                    "whereas the ONT assembly panel uses Flye. The HiFi chain "
                    "reaches the same read-generation and evaluation stages as "
                    "the ONT chain.",
                ],
            ),
            (
                "Quality-value caveat",
                [
                    "The public E. coli HiFi FASTQ used for the cross-platform "
                    "check has an encoding-limited quality distribution: raw "
                    "ASCII quality mean 113.28, Phred-33 mean 80.28 and pysam "
                    "mean 80.08. This is not an off-by-33 script error, but "
                    "absolute QV agreement from this source is reported as "
                    "indicative rather than as a typical benchmark.",
                ],
            ),
            (
                "Human-HiFi status",
                [
                    "The optional human-HiFi source ERR13110527 was incomplete "
                    "and was archived rather than used. Any future human-HiFi "
                    "branch must re-download and validate that source before "
                    "running the analysis.",
                ],
            ),
        ],
        table=[
            ["HiFi measurement", "Value"],
            ["Raw ASCII quality mean", "113.28"],
            ["Phred-33 mean", "80.28"],
            ["pysam mean", "80.08"],
        ],
    )

    print(f"ADDITIONAL_FILES={output}")


if __name__ == "__main__":
    main()
