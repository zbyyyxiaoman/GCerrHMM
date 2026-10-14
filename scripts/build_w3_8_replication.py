#!/usr/bin/env python3
"""Build W3.8 by adding the second A. thaliana replication cohort."""

from __future__ import annotations

import argparse
import json
from copy import deepcopy
from pathlib import Path

from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt
from docx.text.paragraph import Paragraph


def set_text(paragraph: Paragraph, text: str) -> Paragraph:
    if paragraph.runs:
        paragraph.runs[0].text = text
        for run in paragraph.runs[1:]:
            run.text = ""
    else:
        paragraph.add_run(text)
    return paragraph


def find_paragraph(document: Document, prefix: str) -> Paragraph:
    matches = [p for p in document.paragraphs if p.text.strip().startswith(prefix)]
    if len(matches) != 1:
        raise RuntimeError(f"expected one paragraph starting with {prefix!r}, got {len(matches)}")
    return matches[0]


def insert_after(anchor: Paragraph, text: str, template: Paragraph) -> Paragraph:
    new_element = deepcopy(template._p)
    anchor._p.addnext(new_element)
    new_paragraph = Paragraph(new_element, anchor._parent)
    set_text(new_paragraph, text)
    return new_paragraph


def set_cell_text(cell, text: str, *, bold: bool = False) -> None:
    cell.text = ""
    paragraph = cell.paragraphs[0]
    paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
    paragraph.paragraph_format.space_before = Pt(0)
    paragraph.paragraph_format.space_after = Pt(0)
    run = paragraph.add_run(text)
    run.bold = bold
    run.font.name = "Times New Roman"
    run.font.size = Pt(8.5)


def set_table_geometry(table, widths: list[int]) -> None:
    total = sum(widths)
    tbl_pr = table._tbl.tblPr
    for tag in ("tblW", "tblInd"):
        old = tbl_pr.find(qn(f"w:{tag}"))
        if old is not None:
            tbl_pr.remove(old)
    width = OxmlElement("w:tblW")
    width.set(qn("w:w"), str(total))
    width.set(qn("w:type"), "dxa")
    tbl_pr.append(width)
    indent = OxmlElement("w:tblInd")
    indent.set(qn("w:w"), "120")
    indent.set(qn("w:type"), "dxa")
    tbl_pr.append(indent)
    grid = table._tbl.find(qn("w:tblGrid"))
    if grid is not None:
        table._tbl.remove(grid)
    grid = OxmlElement("w:tblGrid")
    for value in widths:
        col = OxmlElement("w:gridCol")
        col.set(qn("w:w"), str(value))
        grid.append(col)
    table._tbl.insert(1, grid)
    for row in table.rows:
        for cell, value in zip(row.cells, widths):
            tc_pr = cell._tc.get_or_add_tcPr()
            old = tc_pr.find(qn("w:tcW"))
            if old is not None:
                tc_pr.remove(old)
            tc_width = OxmlElement("w:tcW")
            tc_width.set(qn("w:w"), str(value))
            tc_width.set(qn("w:type"), "dxa")
            tc_pr.append(tc_width)


def format_p(value: float) -> str:
    value = float(value)
    return "<0.0001" if value < 0.0001 else f"{value:.4f}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--stats", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    payload = json.loads(args.stats.read_text())
    second = payload["second_cohort"]
    first = payload["first_cohort"]
    pooled = payload["all_seed_sensitivity"]

    mean = float(second["mean_delta"])
    low = float(second["ci95_low"])
    high = float(second["ci95_high"])
    p_value = float(second["paired_t_p"])
    pooled_mean = float(pooled["mean_delta"])
    pooled_low = float(pooled["ci95_low"])
    pooled_high = float(pooled["ci95_high"])
    pooled_p = float(pooled["paired_t_p"])
    supported = bool(second["supported"])

    document = Document(args.input)
    body_template = find_paragraph(document, "GCerrHMM is an error-state HMM")

    signpost_replacements = {
        "Results are organised by the claim each block can support.":
        "The results follow the order in which the evidence addresses the "
        "central question.",
        "This block defines the model and the property the paper is allowed to "
        "claim.": "We begin with the model, its parameters and the property "
        "that can be measured, without comparing tools.",
        "The primary mechanistic test is internal.": "The key comparison is "
        "controlled.",
        "A directional difference is interpretable only against the precision "
        "of the instrument.": "The measured difference must be interpreted "
        "against the precision of the instrument.",
        "This ablation is independent of the GC claim.": "The GC comparison is "
        "not affected by this ablation.",
        "The previous blocks test the model mechanism.": "The remaining "
        "results address a different question.",
        "This block separates read-sampling uncertainty from between-window "
        "structural uncertainty":
        "We separate read-sampling uncertainty from between-window structural "
        "uncertainty",
        "The remaining results answer a different question. This block asks a "
        "comparative question:":
        "The remaining results address a comparative question:",
        "The remaining results address a different question.": "",
        "This block asks a comparative question:":
        "We next ask a comparative question:",
        "so this block supports directional scope":
        "so this analysis supports directional scope",
        "uncertainty and states the depth":
        "uncertainty and state the depth",
        "This block therefore isolates the contribution of GC conditioning":
        "The design therefore isolates the contribution of GC conditioning",
        "whatshap": "WhatsHap",
        "GC-stratified error instrument, and more rigorous long-read benchmark "
        "references [10, 11].":
        "GC-stratified error instrument. More rigorous long-read benchmark "
        "references are also needed [10, 11].",
    }
    for paragraph in document.paragraphs:
        for old, new in signpost_replacements.items():
            if old in paragraph.text:
                set_text(paragraph, paragraph.text.replace(old, new, 1))
    for paragraph in document.paragraphs:
        if paragraph.text.strip() == "Summary of evidence and claim scope":
            set_text(paragraph, "Summary of evidence and scope")
        elif paragraph.text.strip() == (
            "Instrument-claim alignment as a reporting principle"
        ):
            set_text(paragraph, "Matching measurements to the stated claim")

    set_text(
        find_paragraph(
            document,
            "The controlled ablation establishes an effect in one structured genome.",
        ),
        "The controlled ablation establishes an effect in one structured "
        "genome. The next analysis asks whether the direction generalises "
        "across six genomes and whether it satisfies the pre-registered "
        "cross-species gate.",
    )
    set_text(
        find_paragraph(
            document,
            "We begin with the model",
        ),
        "We begin with the model, its parameters and the property that can be "
        "measured, without comparing tools. This section states the transition "
        "model, the GC binning and the controlled contrast tested below.",
    )
    set_text(
        find_paragraph(
            document,
            "The results follow the order",
        ),
        "We first define the model and the fidelity problem, then isolate GC "
        "conditioning with a one-bin control. The next analysis tests whether "
        "the effect generalises across six genomes under the pre-registered "
        "gate. We then quantify measurement uncertainty, report an independent "
        "indel ablation, compare descriptive and downstream performance with "
        "existing simulators, and close with computational scaling and the "
        "boundary of the evidence. The composite score is not used as evidence "
        "for GC conditioning.",
    )

    set_text(
        find_paragraph(document, "Article type:"),
        "Article type: Research article | Target journal: BMC Bioinformatics | "
        'Collection "Simulated and synthetic data in bioinformatics: statistical '
        'methods, models, and applications" | W3.8 major-revision manuscript, '
        "2026-10-14.",
    )

    abstract_result = (
        "confirmed the single-species effect"
        if supported
        else "did not confirm the single-species effect"
    )
    set_text(
        find_paragraph(document, "We present GCerrHMM"),
        "We present GCerrHMM, a transition-based error HMM trained on real ONT "
        "alignments whose state-transition model is conditioned on local GC "
        "content, with a reproducible fidelity harness. The model is evaluated "
        "with a controlled one-bin versus GC-aware ablation, a six-genome "
        "directional panel, and coverage-matched human chr21 downstream "
        "layers. GC-aware simulation moved the local GC-error curve in the "
        "expected direction for four of six genomes in both replicates; the "
        "strict cross-species gate was not met, so the cross-species result "
        "remains directional. In A. thaliana, five or more bins raised curve "
        "correlation from r = 0.779 (SD 0.062) to r >= 0.960. A second, "
        f"pre-specified six-seed cohort {abstract_result} "
        f"(mean delta r = {mean:.4f}; 95% CI {low:.4f} to {high:.4f}; "
        f"paired t-test p = {p_value:.4f}), while the original cross-species "
        "gate remains unchanged. The harness separates read-sampling from "
        "between-window uncertainty, and process-level scaling preserves "
        "FASTQ identity while reducing wall time.",
    )
    conclusion = (
        "A pre-specified single-species cohort confirms the A. thaliana "
        "GC-conditioning effect"
        if supported
        else "The second pre-specified single-species cohort did not confirm "
        "the A. thaliana GC-conditioning effect"
    )
    set_text(
        find_paragraph(document, "GCerrHMM provides a trainable"),
        "GCerrHMM provides a trainable GC-aware error model, a reproducible "
        "evaluation harness, and a process-scalable simulator. "
        f"{conclusion}, while cross-species generalization remains directional. "
        "The study shows that simulator benchmarking must align the measurement "
        "instrument with the claim and report both statistical and "
        "computational performance.",
    )

    confirmation = find_paragraph(document, "Confirmatory replication.")
    phase2_result = (
        "confirmed the single-species effect"
        if supported
        else "did not confirm the single-species effect"
    )
    phase2_clause = "was confirmed" if supported else "was not confirmed"
    phase2_text = (
        "Second pre-specified replication. A second independent A. thaliana "
        "cohort was fixed before running with six new seeds (77, 88, 99, 111, "
        "122 and 133) and the same primary endpoint, the paired difference in "
        f"GC-curve Pearson r between 5-bin and 1-bin models. The cohort "
        f"{phase2_result} (mean delta r = {mean:.4f}; "
        f"95% CI {low:.4f} to {high:.4f}; paired t-test p = {p_value:.4f}). "
        "The first independent cohort (seeds 44, 55 and 66) remains reported "
        f"separately as directionally positive but not confirmatory "
        f"(p = {float(first['paired_t_p']):.4f}). Across all 12 seeds from the "
        f"discovery, first-replication and second-replication cohorts, the "
        f"pooled sensitivity analysis gave mean delta r = {pooled_mean:.4f} "
        f"(95% CI {pooled_low:.4f} to {pooled_high:.4f}; "
        f"p {format_p(pooled_p)}). The original cross-species gate is not reopened "
        "or retroactively changed. These cohorts are summarized in Table 6.",
    )
    insert_after(confirmation, phase2_text, body_template)

    table_caption = insert_after(
        find_paragraph(document, "Second pre-specified replication."),
        "Table 6. Pre-specified A. thaliana replication cohorts.",
        body_template,
    )

    # Move the old depth caption to Table 8; Tables 7 and 8 keep their
    # narrative slots (six-species panel and depth), while the resource and
    # relationship tables move to 9 and 10.
    depth_caption = find_paragraph(document, "Table 6. Depth dependence")
    set_text(depth_caption, "Table 8. Depth dependence of the read-sampling component.")
    for caption in document.paragraphs:
        if caption.text.strip().startswith("Table 8. Computational resource"):
            set_text(
                caption,
                "Table 9. Computational resource benchmark on the E. coli "
                "reference and real ONT alignment.",
            )
            break
    for caption in document.paragraphs:
        if caption.text.strip().startswith("Table 9. Relationship"):
            set_text(
                caption,
                "Table 10. Relationship to representative long-read simulators.",
            )
            break
    for paragraph in document.paragraphs:
        if paragraph.text.strip().startswith("Table 8 reports"):
            set_text(
                paragraph,
                paragraph.text.replace("Table 8 reports", "Table 9 reports", 1),
            )
            break

    table = document.add_table(rows=5, cols=6)
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = False
    headers = ["Cohort", "Seeds", "n", "Mean delta r", "95% CI", "p"]
    for index, header in enumerate(headers):
        set_cell_text(table.rows[0].cells[index], header, bold=True)
    rows = [
        ("Discovery", "11, 22, 33", "3",
         f"{float(payload['discovery_cohort']['mean_delta']):.4f}",
         f"[{float(payload['discovery_cohort']['ci95_low']):.4f}, "
         f"{float(payload['discovery_cohort']['ci95_high']):.4f}]",
         f"{float(payload['discovery_cohort']['paired_t_p']):.4f}"),
        ("First replication", "44, 55, 66", "3",
         f"{float(first['mean_delta']):.4f}",
         f"[{float(first['ci95_low']):.4f}, {float(first['ci95_high']):.4f}]",
         f"{float(first['paired_t_p']):.4f}"),
        ("Second replication", "77, 88, 99, 111, 122, 133", "6",
         f"{mean:.4f}", f"[{low:.4f}, {high:.4f}]", f"{p_value:.4f}"),
        ("Pooled sensitivity", "11, 22, 33, 44, 55, 66, 77, 88, 99, 111, 122, 133", "12",
         f"{pooled_mean:.4f}", f"[{pooled_low:.4f}, {pooled_high:.4f}]",
         format_p(pooled_p)),
    ]
    for row_index, values in enumerate(rows, start=1):
        for column_index, value in enumerate(values):
            set_cell_text(table.rows[row_index].cells[column_index], value)
    set_table_geometry(table, [1000, 3200, 450, 1100, 2400, 1210])
    table_caption._p.addnext(table._tbl)

    discussion = find_paragraph(document, "The six-genome panel spans")
    set_text(
        discussion,
        "The six-genome panel spans GC-uniform and strongly GC-heterogeneous "
        "genomes and supports a boundary hypothesis: the strongest positive "
        "case was H. sapiens chr21; A. thaliana was mixed across the discovery "
        f"route panel, and the 1-bin versus 5-bin contrast {phase2_clause} in a "
        "second pre-specified six-seed cohort; D. melanogaster remained mixed. "
        "The pre-registered directional gate still defines the scope of the "
        "cross-species claim, and the human panel remains a two-replicate "
        "supporting analysis. A future pre-registered panel with more genomes "
        "and independent read sets should test whether the model-implied "
        "GC-error amplitude predicts when GC-aware training is beneficial.",
    )

    limitations = find_paragraph(document, "First, the cross-species GC-conditional")
    set_text(
        limitations,
        "First, the cross-species GC-conditional effect remains directional "
        "evidence with effect sizes across 10-30x. The second pre-specified "
        f"A. thaliana 1-bin versus 5-bin cohort {phase2_clause}, while the "
        "first independent cohort was directionally "
        "positive but not confirmatory and the original cross-species gate "
        "remains unmet. The measured depth budget is reported in the Abstract. "
        + "Second," + limitations.text.split("Second,", 1)[1],
    )

    simulation = find_paragraph(
        document,
        "The external cross-tool panel uses one simulated replicate",
    )
    set_text(
        simulation,
        simulation.text.replace(
            "and a separate prospectively specified A. thaliana 1-bin versus "
            "5-bin replication uses seeds 44, 55 and 66 at the same coverage.",
            "a separate prospectively specified A. thaliana 1-bin versus "
            "5-bin replication uses seeds 44, 55 and 66 at the same coverage; "
            "a second pre-specified six-seed cohort uses seeds 77, 88, 99, "
            "111, 122 and 133 with the same endpoint and one-time analysis.",
        ),
    )

    ablation_design = find_paragraph(
        document,
        "Two ablation axes were pre-registered.",
    )
    set_text(
        ablation_design,
        "Two ablation axes were pre-registered. GC-bin granularity used K in "
        "{1, 5, 10, 20} with profile-matched training and three seeds per "
        "setting, evaluated on both the composite and the GC-stratified "
        "instrument; the pre-registered gate defines the scope of the "
        "cross-species significance claim. Training coverage was assessed on "
        "matched panels across a coverage grid with the composite. The "
        "observed effects are reported with effect sizes and the measurement "
        "budget. A separate pre-specified confirmation program tested only "
        "the A. thaliana 1-bin versus 5-bin curve-correlation contrast. The "
        "first cohort used seeds 44, 55 and 66. Because it was directionally "
        "positive but did not meet the p < 0.05 threshold, a second "
        "independent cohort using seeds 77, 88, 99, 111, 122 and 133 was "
        "fixed before running, with the same endpoint and a one-time analysis. "
        "Neither cohort replaces or reopens the original cross-species "
        "decision rule.",
    )

    set_text(
        find_paragraph(document, "Additional file 3."),
        "Additional file 3. File name: Additional_file_3_gc_bin_ablation.xlsx. "
        "File format: XLSX. Title: GC-bin ablation and A. thaliana replication "
        "full table. Description: curve r and MAD by species, bin setting and "
        "seed, including the first and second independent A. thaliana cohorts "
        "and the pooled 12-seed sensitivity analysis.",
    )

    for section in document.sections:
        for header_paragraph in section.header.paragraphs:
            if "W3.7" in header_paragraph.text:
                set_text(
                    header_paragraph,
                    header_paragraph.text.replace("W3.7", "W3.8"),
                )
    document.core_properties.comments = (
        "Local W3.8 revision adding the second pre-specified A. thaliana "
        "replication cohort."
    )
    document.core_properties.last_modified_by = "GCerrHMM project"

    args.output.parent.mkdir(parents=True, exist_ok=True)
    document.save(args.output)
    print(args.output)


if __name__ == "__main__":
    main()
