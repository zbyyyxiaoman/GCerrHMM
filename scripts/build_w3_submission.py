#!/usr/bin/env python3
"""Convert the W2 layout draft into a BMC submission-format manuscript."""

from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

from docx import Document
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Pt
from docx.table import Table
from docx.text.paragraph import Paragraph


FIGURE_MAP = {
    "background_figure_b1_gc_context": "Figure1_gc_context",
    "background_figure_b2_simulator_landscape": "Figure2_simulator_landscape",
    "background_figure_b3_context_layers": "Figure3_context_layers",
    "figure1_overview_innovation": "Figure4_overview_innovation",
    "figure2_reads_cross_tool": "Figure5_reads_cross_tool",
    "figure3_alignment": "Figure6_alignment",
    "figure4_variant_calling": "Figure7_variant_calling",
    "figure5_assembly_phasing": "Figure8_assembly_phasing",
    "figure6_ablation": "Figure9_ablation",
}


def find_paragraph(document: Document, needle: str) -> Paragraph:
    for paragraph in document.paragraphs:
        if needle.lower() in paragraph.text.lower():
            return paragraph
    raise RuntimeError(f"paragraph not found: {needle}")


def set_paragraph_text(paragraph: Paragraph, text: str) -> None:
    """Replace paragraph text while retaining the paragraph object."""
    if paragraph.runs:
        paragraph.runs[0].text = text
        for run in paragraph.runs[1:]:
            run.text = ""
    else:
        paragraph.text = text


def create_paragraph_after(
    anchor: Paragraph,
    text: str = "",
    style: str | None = None,
) -> Paragraph:
    element = OxmlElement("w:p")
    anchor._p.addnext(element)
    paragraph = Paragraph(element, anchor._parent)
    paragraph.text = text
    if style:
        try:
            paragraph.style = style
        except KeyError:
            pass
    return paragraph


def renumber_figure_text(text: str) -> str:
    """Renumber Background B1-B3 and main Figure 1-6 into one sequence."""
    placeholders: dict[str, str] = {}

    def background_replacement(match: re.Match[str]) -> str:
        key = f"__BG_FIGURE_{match.group(1)}__"
        placeholders[key] = f"Figure {int(match.group(1))}"
        return key

    text = re.sub(
        r"Background Figure B([123])",
        background_replacement,
        text,
    )
    for old_number in range(6, 0, -1):
        new_number = old_number + 3
        ref_key = f"__FIG_REF_{new_number}__"
        figure_key = f"__FIG_NUMBER_{new_number}__"
        placeholders[ref_key] = f"Fig. {new_number}"
        placeholders[figure_key] = f"Figure {new_number}"
        text = re.sub(
            rf"\bFig\.\s*{old_number}(?!\d)",
            ref_key,
            text,
        )
        text = re.sub(
            rf"\bFigure {old_number}(?!\d)",
            figure_key,
            text,
        )
    for key, value in placeholders.items():
        text = text.replace(key, value)
    return text


def remove_embedded_figures(document: Document) -> dict[str, str]:
    """Remove review-only images and collect Background captions."""
    paragraphs = list(document.paragraphs)
    remove_ids: set[int] = set()
    background_captions: dict[str, str] = {}
    for index, paragraph in enumerate(paragraphs):
        if "w:drawing" not in paragraph._p.xml:
            continue
        remove_ids.add(id(paragraph._p))
        if index + 1 >= len(paragraphs):
            continue
        following = paragraphs[index + 1]
        following_text = following.text.strip()
        match = re.match(
            r"^Background Figure (B[123])\.\s*(.+)$",
            following_text,
            flags=re.DOTALL,
        )
        if match:
            background_captions[match.group(1)] = re.sub(
                r"\s+", " ", match.group(2)
            ).strip()
        if re.match(r"^(?:Background )?Figure ", following_text):
            remove_ids.add(id(following._p))
    for paragraph in paragraphs:
        if id(paragraph._p) in remove_ids:
            paragraph._p.getparent().remove(paragraph._p)
    return background_captions


def append_to_paragraph(
    document: Document,
    prefix: str,
    sentence: str,
) -> None:
    for paragraph in document.paragraphs:
        if paragraph.text.strip().startswith(prefix):
            set_paragraph_text(
                paragraph,
                paragraph.text.rstrip() + " " + sentence.strip(),
            )
            return
    raise RuntimeError(f"paragraph not found for citation: {prefix}")


def repair_incomplete_sentences(document: Document) -> None:
    """Repair line-break artefacts carried over from the W1 skeleton."""
    paragraphs = list(document.paragraphs)
    for index, paragraph in enumerate(paragraphs):
        text = paragraph.text.strip()
        if text.startswith("We present GCerrHMM") and text.endswith("meeting"):
            set_paragraph_text(
                paragraph,
                text[:-len("meeting")].rstrip()
                + " meeting the pre-registered top-2 criterion in both "
                "anchored layers.",
            )
        elif text.startswith("1.01 at 10x"):
            if text.endswith("meeting"):
                set_paragraph_text(
                    paragraph,
                    text[:-len("meeting")].rstrip()
                    + " meeting the pre-registered top-2 criterion in "
                    "both anchored layers.",
                )
            if index > 0 and paragraphs[index - 1].text.strip().startswith(
                "We present GCerrHMM"
            ):
                set_paragraph_text(
                    paragraphs[index - 1],
                    paragraphs[index - 1].text.rstrip() + " "
                    + paragraph.text.strip(),
                )
                paragraph._p.getparent().remove(paragraph._p)
        elif text.startswith(
            "On the coverage-matched chr21 30x delta-to-real ranking"
        ):
            if not text.endswith((".", "!", "?")):
                set_paragraph_text(paragraph, text + ".")
        elif text.startswith("Two downstream layers carry"):
            set_paragraph_text(
                paragraph,
                "Two downstream layers carry a real-data anchor measured "
                "through the identical pipeline: alignment and assembly. "
                "On the assembly panel (Flye for ONT, hifiasm for HiFi; "
                "PBSim3-errhmm omitted because no assembly run exists), "
                "assembly metrics were compared against the real 30x ONT "
                "anchor. The assembly and phasing summaries are shown in "
                "Fig. 8.",
            )


def replace_prefix(
    document: Document,
    prefix: str,
    replacement: str,
) -> None:
    for paragraph in document.paragraphs:
        if paragraph.text.strip().startswith(prefix):
            set_paragraph_text(paragraph, replacement)
            return
    raise RuntimeError(f"paragraph not found for replacement: {prefix}")


def revise_scientific_text(document: Document) -> None:
    """Apply W3.1 corrections justified by code, data and pre-registration."""
    replace_prefix(
        document,
        "With a single insertion state",
        (
            "With a single insertion state, a one-base and a fifty-base "
            "insertion collapse to the same state, the I-to-I transition is "
            "never observed, and the learned model emits every insertion as "
            "exactly 1 bp. On the human chr21 training alignment the mean "
            "insertion length is 3.83 bp, and the single-state encoding "
            "produces 0.04% inserted bases against 3.78% under the same "
            "CIGAR-derived M-base denominator. The run-length encoding "
            "(I1-I3, I4+; D1-D3, D4+, with a heavy-tailed length distribution "
            "for the 4+ states) moves simulated chr21 reads from 0.04% to "
            "2.54% inserted bases and from a fixed 1.00 bp to 2.83 bp mean "
            "insertion length, and raises the total error rate from 6.5% to "
            "8.8% against a real 14.5%. As a single-variable ablation on "
            "E. coli, the encoding improves the reads-level composite from "
            "83.5 to 87.6, driven by a 16-point gain in the k-mer sub-score. "
            "The residual gap to real data is concentrated in deletion run "
            "length and in read segments that alignment-based training "
            "cannot observe."
        ),
    )
    replace_prefix(
        document,
        "Indel run length is carried",
        (
            "Indel run length is carried in the state name, for insertions "
            "and deletions alike. A run of one, two or three bases enters "
            "I1/I2/I3 or D1/D2/D3; longer runs share I4+/D4+, whose length "
            "distribution is learned separately and is heavy-tailed. Among "
            "chr21 insertion runs of 4 bp or more, 35% are exactly 4 bp, 52% "
            "are 5-9 bp, 10% are 10-49 bp and 3.5% are at least 50 bp, with "
            "a mean of about 241 bp in that final tail. With a single I state "
            "the I-to-I transition is never observed and the model emits "
            "every insertion as exactly 1 bp; the quantitative consequences "
            "are reported in Results."
        ),
    )
    replace_prefix(
        document,
        "We present GCerrHMM",
        (
            "We present GCerrHMM, a transition-based error HMM trained on "
            "real ONT alignments whose state-transition model is conditioned "
            "on local GC content, and a reproducible fidelity harness. "
            "Across six genomes the GC-aware "
            "model was directionally closer to the real GC-error curve than "
            "its one-bin control in both replicates for four species, "
            "including the GC-heterogeneous H. sapiens chr21; "
            "A. thaliana and D. melanogaster were mixed. We report effect "
            "sizes and the precision of the between-window comparison at "
            "10x, 20x and 30x. In a "
            "separate depth grid the read-sampling SD ranged from 0.009 to "
            "0.194 across 5-15x (6/3/2 subsets); the 15x estimate is "
            "suggestive rather than precise. The composite distributional "
            "score, which contains no GC-conditional term, is reported as a "
            "descriptive summary; the GC-stratified instrument carries the "
            "GC-specific comparison. On "
            "the coverage-matched chr21 30x panel, GCerrHMM ranked second on "
            "alignment identity and first on assembly Borda distance, "
            "meeting the pre-registered top-2 criterion in both anchored "
            "layers."
        ),
    )
    replace_prefix(
        document,
        "GCerrHMM is a trainable, GC-aware error",
        (
            "GCerrHMM provides a trainable GC-aware error model and a "
            "reproducible evaluation harness. Its GC-conditional effect is "
            "measurable in structure; we report the effect size, the separate "
            "sampling and structural uncertainties, and the depth needed for "
            "a definitive cross-species test. The study "
            "also shows that simulator benchmarking requires "
            "instrument-claim alignment: a fidelity score can test only the "
            "property it contains."
        ),
    )
    replace_prefix(
        document,
        "GCerrHMM is a trainable, GC-aware error HMM",
        (
            "GCerrHMM is a trainable GC-aware error model for long-read "
            "simulation, with a length-aware indel encoding, competitive "
            "compositional fidelity, and downstream consistency that ranks "
            "first or second among tested simulators on real-anchored "
            "layers. Its GC-conditional effect is measurable in structure; "
            "we report effect sizes, the separate sampling and structural "
            "uncertainties, and the depth budget required for a definitive "
            "cross-species test. The study demonstrates that simulator "
            "benchmarking requires instrument-claim alignment: a fidelity "
            "score can test only the property it contains. The fail-loud "
            "harness, instrument audit and pre-registered decision framework "
            "are released with the model."
        ),
    )
    replace_prefix(
        document,
        "The six-genome panel was not chosen arbitrarily",
        (
            "The six-genome panel spans GC-uniform and strongly "
            "GC-heterogeneous genomes and provides directional support for a "
            "boundary hypothesis: the strongest positive case was "
            "H. sapiens chr21, while A. thaliana and D. melanogaster were "
            "mixed. The pre-registered directional gate defines the scope of "
            "the cross-species claim, and the human panel remains a "
            "two-replicate supporting analysis. A future pre-registered panel with more "
            "genomes and independent read sets should test whether the "
            "model-implied GC-error amplitude predicts when GC-aware "
            "training is beneficial."
        ),
    )
    replace_prefix(
        document,
        "Across the six-species route panel",
        (
            "Across the six-species route panel, the k-mer statistic is not "
            "a GC-stratified instrument: 11 of 12 route/replicate "
            "comparisons had a difference interval containing zero, and the "
            "single interval excluding zero favoured the 1-bin control "
            "without reproducing. That k-mer panel is therefore reported as "
            "not discriminable for distributional composition. The "
            "GC-stratified curve comparison is reported separately in "
            "Table 7: the GC-aware model is directionally favoured in both "
            "replicates for four of six genomes, while A. thaliana and "
            "D. melanogaster are mixed. The pre-registered double-species "
            "gate defines the scope of the cross-species claim; we report "
            "the observed directional evidence and effect sizes."
        ),
    )
    replace_prefix(
        document,
        "On the six-species panel",
        (
            "On the six-species panel, the GC-aware model exceeded its "
            "one-bin control in both replicates for E. coli "
            "(0.458/0.589 vs. -0.738/0.505), S. cerevisiae "
            "(0.471/0.535 vs. 0.230/0.174), M. musculus chr19 "
            "(0.922/0.794 vs. 0.460/0.650) and H. sapiens chr21 "
            "(0.917/0.861 vs. -0.579/-0.091). A. thaliana and "
            "D. melanogaster were mixed. The strongest directional signal "
            "came from H. sapiens chr21; the pre-registered double-species "
            "gate defines the scope of the cross-species claim, and the "
            "human panel remains a two-replicate supporting analysis. In "
            "the 10x cross-tool panel, the comparison provides a "
            "conservative reference point; the GC-stratified six-species "
            "panel is reported in Table 7."
        ),
    )
    replace_prefix(
        document,
        "For phasing, real HG002 HiFi reads",
        (
            "For phasing, real HG002 HiFi reads processed with WhatsHap "
            "achieved a phased fraction of 0.9655 across 179 blocks with "
            "block N50 of 177,187 bp. Switch error is not reported because "
            "no phased truth set is bundled with the panel. The spike-in "
            "donor designs are haploid and carry no heterozygous sites, so "
            "no simulated-read phasing arm is reported. The phasing panel is "
            "therefore a real-data quality check, not a cross-tool "
            "comparison."
        ),
    )
    replace_prefix(
        document,
        "Panel models used the default minimum MAPQ",
        (
            "Panel models used the default minimum MAPQ of 20, except for "
            "the human chr21 30x model. Chr21 was the only species whose "
            "read pool fell below the 20,000-read cap: 58% of its 8,189 "
            "primary reads mapped below MAPQ 20, so applying the default "
            "retained 3,438 reads and biased the learned error spectrum "
            "(event rate 0.0459 at MAPQ 20 vs 0.0713 at MAPQ 0). The chr21 "
            "model therefore uses minimum MAPQ 0; the other species retained "
            "MAPQ 20 because their larger pools were capped at 20,000 reads "
            "and the threshold changed the retained event rate by less than "
            "0.1 percentage points. This is a platform-specific correction "
            "to an alignment-depth artefact, not a different model class. "
            "Hierarchical smoothing used prior strength 10.0 with "
            "pseudocount 1.0."
        ),
    )
    replace_prefix(
        document,
        "The species table is generated",
        (
            "The species table is generated from "
            "config/species_config_v2.json and config/data_sources.json; "
            "Fig. 4c lists reference and ONT run accessions. Training uses "
            "the real alignment panels listed there, whereas evaluation uses "
            "matched simulated FASTQ at the coverage stated in each Results "
            "table; no simulated reads are used to train the reported model."
        ),
    )
    replace_prefix(
        document,
        "The GC claim is tested",
        (
            "The GC claim is tested with an aligned instrument: per-window "
            "GC bins (100 bp, equal width, the trainer's convention) and the "
            "error rate of aligned bases within each bin from the CIGAR with "
            "soft-clipped bases excluded. Real and simulated sets are aligned "
            "with the same mapper; the simulated per-bin error profile is "
            "compared with the real one by Pearson and Spearman correlation "
            "and mean absolute deviation. A one-bin model cannot generate a "
            "GC-stratified error profile, so its finite-sample correlation "
            "across GC bins reflects alignment between the global error "
            "level and the real curve, not GC-conditional fidelity; the "
            "pre-registered decision compares the one-bin control against "
            "the GC-aware model rather than testing whether the one-bin "
            "correlation is identically zero. The k-mer term correlates "
            "frequency vectors over the intersection of per-side top-10,000 "
            "k-mer maps (k = 21, Pearson); the sample was increased from 250 "
            "to 1000 reads per side with a paired bootstrap over reads after "
            "the 250-read estimator showed a replicate spread of up to 16.7 "
            "points."
        ),
    )
    replace_prefix(
        document,
        "The composite aggregates read-length",
        (
            "The composite aggregates read-length, QV, GC-content and "
            "k-mer terms; it is designed for marginal distributional "
            "similarity and is reported as a descriptive summary. The GC "
            "claim is evaluated on a GC-stratified error curve - "
            "per-window GC bins (100 bp, equal width) and the error rate of "
            "aligned bases within each bin - compared between real and "
            "simulated reads by Pearson correlation and mean absolute "
            "deviation, with explicit intervals (Methods). A one-bin model "
            "cannot generate a GC-stratified profile; its finite-sample "
            "correlation therefore reflects alignment between the global "
            "error level and the real curve, not GC-conditional capability. "
            "The pre-registered decision compares the one-bin control "
            "against the GC-aware model across bins and species rather than "
            "testing whether a near-zero one-bin correlation is a separate "
            "hypothesis."
        ),
    )
    replace_prefix(
        document,
        "The 95% between-window bootstrap interval",
        (
            "The 95% between-window bootstrap interval on the GC-curve "
            "correlation was 1.01 at 10x, 0.96 at 20x and 0.96 at 30x. It "
            "does not shrink with sequencing depth and is a structural "
            "uncertainty of the window-level curve, not a power statement. "
            "In the separate depth grid, the read-sampling SD was 0.176 at "
            "5x (6 disjoint subsets), 0.194 at 10x (3 subsets) and 0.009 at "
            "15x (2 halves). The 15x estimate is therefore suggestive rather "
            "than precise; it indicates that an effect larger than about "
            "0.05 correlation units would have been detectable above 15x, "
            "while the between-window interval remains wide."
        ),
    )
    for paragraph in document.paragraphs:
        if paragraph.text.strip() == (
            "A boundary condition: the benefit scales with GC-structure "
            "heterogeneity"
        ):
            set_paragraph_text(paragraph, "A testable boundary hypothesis")
        if paragraph.text.strip() == "Relation to existing long-read simulators":
            try:
                paragraph.style = "Heading 2"
            except KeyError:
                pass
        if paragraph.text.strip() == (
            "Six-species panel: not discriminable, with an interpretable "
            "boundary"
        ):
            set_paragraph_text(
                paragraph,
                "Six-species panel: distributional and GC-stratified results",
            )
        if paragraph.text.strip().startswith("We set out to ask"):
            set_paragraph_text(
                paragraph,
                "We asked whether conditioning a long-read error model on "
                "local GC content improves simulation fidelity and whether "
                "standard evaluation metrics can detect that improvement. "
                "We found measurable GC-conditional error structure and "
                "directionally closer GC-error curves in four of six "
                "genomes; the pre-registered gate defines the scope of the "
                "cross-species claim. The broader lesson is that a fidelity "
                "score can test "
                "only the property it contains; evaluation instruments and "
                "claims must be aligned.",
            )
        if paragraph.text.strip() == (
            "The composite does not separate bin settings; the stratified "
            "instrument gives directional evidence"
        ):
            set_paragraph_text(
                paragraph,
                "Composite separability and stratified evidence",
            )
    replace_prefix(
        document,
        "For each 100-bp reference window",
        (
            "For each 100-bp reference window, local GC content is assigned "
            "to one of K equal-width bins. The model estimates "
            "P(state at i+1 | state at i, GC bin) from the CIGAR path "
            "together with a nucleotide substitution matrix. Conditioning "
            "therefore applies to the state-transition matrix; the "
            "substitution matrix is global rather than GC-specific. The "
            "default state space is M, S, I1, I2, I3, I4+, D1, D2, D3, D4+; "
            "a simplified state space merges the indel states into single I "
            "and D states and serves as an ablation control. The one-bin "
            "control (K = 1) is trained on the same aligned reads with the "
            "same MAPQ threshold, GC window and read budget."
        ),
    )
    replace_prefix(
        document,
        "The model-implied error rate per GC bin",
        (
            "The model-implied error rate per GC bin recovers the "
            "qualitative shape of the real data where such structure "
            "exists: in A. thaliana, whose genome carries strong local GC "
            "heterogeneity, the trained model reproduces a sharp error peak "
            "in the highest-GC bin (19.5% in bin 7), while the curves for "
            "the GC-uniform E. coli and S. cerevisiae genomes are flat "
            "(Fig. 4b). Fig. 4b is model-implied; the real-versus-simulated "
            "GC-error curves used to motivate the model are shown in Fig. 1b. "
            "The six-genome panel and its data sources are summarised in "
            "Fig. 4c."
        ),
    )
    replace_prefix(
        document,
        "First, the GC-conditional advantage",
        (
            "First, the GC-conditional effect is reported as directional "
            "evidence with effect sizes across 10-30x; the pre-registered "
            "gate defines the scope of the cross-species claim, and the "
            "measured depth budget is reported in the Abstract. Second, the "
            "Level-1 composite measures distributional similarity only; "
            "conclusions about GC conditioning rest on the stratified-error "
            "instrument and the controlled bin ablation. Third, the k-mer "
            "sub-score has little discriminating power for genomes with a "
            "small shared k-mer intersection (D. melanogaster r = 0.12-0.35 "
            "vs. E. coli 0.86-0.94) and a replicate spread of up to 16.7 "
            "points. Fourth, the variant-calling layer uses a synthetic "
            "truth set and has no real-data anchor; re-anchoring it to GIAB "
            "HG002 high-confidence calls is future work, as is the 6.62x "
            "assembly condition where real reads do not assemble. Fifth, "
            "phasing switch error could not be reported for lack of a "
            "bundled phased truth. Sixth, the public HiFi E. coli FASTQ "
            "used for the cross-platform check has an encoding-limited "
            "quality distribution, so HiFi QV agreement is indicative only."
        ),
    )
    replace_prefix(
        document,
        "Here we present GCerrHMM",
        (
            "Here, we present GCerrHMM, an error-state hidden Markov model "
            "trained from real long-read alignments whose transition "
            "structure is conditioned on local GC content, with a "
            "length-aware indel encoding that removes a structural bias we "
            "show is present in single-state encodings. We evaluate it with "
            "a two-level harness (distributional fidelity; downstream task "
            "consistency) across six genomes and four simulators, with "
            "pre-registered decision rules for the central GC claim. We "
            "report the observed directional effects, the precision limits "
            "of the current instrument, and the depth budget required for a "
            "definitive cross-species test. Fig. 3 places the present study "
            "within that context and evaluation chain."
        ),
    )
    replace_prefix(
        document,
        "Two ablation axes were pre-registered",
        (
            "Two ablation axes were pre-registered. GC-bin granularity: K "
            "in {1, 5, 10, 20} with profile-matched training and three seeds "
            "per setting, evaluated on both the composite and the "
            "GC-stratified instrument; the pre-registered gate defines the "
            "scope of the cross-species significance claim. Training "
            "coverage: matched panels across a coverage grid, evaluated on "
            "the composite. The observed effects are reported with effect "
            "sizes and the measurement budget."
        ),
    )

    for paragraph in document.paragraphs:
        text = paragraph.text
        text = text.replace("Fig. 4c now lists", "Fig. 4c lists")
        text = text.replace(" --gcode", " --gc-demo")
        text = text.replace("WhatsApp", "WhatsHap")
        text = text.replace("whatsap", "WhatsHap")
        text = text.replace("badead", "badread")
        text = text.replace("bradread", "badread")
        text = text.replace("GCerHMM", "GCerrHMM")
        text = text.replace("Here we present", "Here, we present")
        if text != paragraph.text:
            set_paragraph_text(paragraph, text)


def normalize_table_labels(document: Document) -> None:
    replacements = {
        "Ecoli": "E. coli",
        "Scerevisiae": "S. cerevisiae",
        "Athaliana": "A. thaliana",
        "Dmelanogaster": "D. melanogaster",
        "Mmusculus_chr19": "M. musculus chr19",
        "Hsapiens_chr21": "H. sapiens chr21",
        "flye": "Flye",
        "GC-conditioned error HMM": "GC-aware error HMM",
    }
    for table in document.tables:
        for row in table.rows:
            for cell in row.cells:
                for paragraph in cell.paragraphs:
                    text = paragraph.text
                    for old, new in replacements.items():
                        text = text.replace(old, new)
                    if text != paragraph.text:
                        set_paragraph_text(paragraph, text)


def find_table_by_caption(
    document: Document,
    caption_prefix: str,
) -> Table:
    for child in document.element.body.iterchildren():
        if not child.tag.endswith("}tbl"):
            continue
        previous = child.getprevious()
        while previous is not None and not previous.tag.endswith("}p"):
            previous = previous.getprevious()
        if previous is None:
            continue
        caption = Paragraph(previous, document).text.strip()
        if caption.startswith(caption_prefix):
            return Table(child, document)
    raise RuntimeError(f"table not found: {caption_prefix}")


def add_note_after_table(
    document: Document,
    caption_prefix: str,
    note: str,
) -> None:
    for child in document.element.body.iterchildren():
        if not child.tag.endswith("}tbl"):
            continue
        previous = child.getprevious()
        while previous is not None and not previous.tag.endswith("}p"):
            previous = previous.getprevious()
        if previous is None:
            continue
        caption = Paragraph(previous, document).text.strip()
        if not caption.startswith(caption_prefix):
            continue
        element = OxmlElement("w:p")
        child.addnext(element)
        paragraph = Paragraph(element, document)
        paragraph.text = note
        try:
            paragraph.style = "Caption"
        except KeyError:
            pass
        return
    raise RuntimeError(f"table not found for note: {caption_prefix}")


def add_real_anchor_row(document: Document) -> None:
    table = find_table_by_caption(document, "Table 4.")
    if any(
        "Real 30x ONT anchor" in cell.text
        for row in table.rows
        for cell in row.cells
    ):
        return
    cells = table.add_row().cells
    values = [
        "H. sapiens chr21",
        "Real 30x ONT anchor",
        "Flye",
        "48",
        "33.50",
        "79.77",
        "1",
    ]
    for cell, value in zip(cells, values):
        cell.text = value


def move_relation_section(document: Document) -> None:
    relation = find_paragraph(
        document,
        "Relation to existing long-read simulators",
    )
    relation_properties = relation._p.get_or_add_pPr()
    relation_style = relation_properties.find(qn("w:pStyle"))
    if relation_style is None:
        relation_style = OxmlElement("w:pStyle")
        relation_properties.insert(0, relation_style)
    relation_style.set(qn("w:val"), "Heading2")
    limitations = find_paragraph(document, "Limitations")
    first_limitation = next(
        (
            paragraph
            for paragraph in document.paragraphs
            if paragraph.text.strip().startswith(
                (
                    "First, the GC-conditional advantage",
                    "First, the GC-conditional effect",
                )
            )
        ),
        None,
    )
    if first_limitation is None:
        raise RuntimeError("first limitations paragraph not found")
    blocks = []
    element = relation._p
    while element is not None and element is not first_limitation._p:
        blocks.append(element)
        element = element.getnext()
    for block in blocks:
        limitations._p.addprevious(block)


def add_table_notes(document: Document) -> None:
    for paragraph in document.paragraphs:
        if paragraph.text.strip().startswith("Table 4."):
            set_paragraph_text(
                paragraph,
                "Table 4. Assembly N50 and reference identity for the "
                "simulated ONT panel and the real 30x ONT anchor "
                "(mean; n=1 per row).",
            )
            break
    add_real_anchor_row(document)
    add_note_after_table(
        document,
        "Table 2.",
        "n = one matched run per tool-species row. PBSim3-errhmm uses the "
        "same errhmm mode across species; its similar identity values are "
        "measured, not imputed, and each row comes from an independent run.",
    )
    add_note_after_table(
        document,
        "Table 3.",
        "R3 is the absolute difference in base identity from the real 30x "
        "ONT anchor. R5 Borda is the mean of the rank by absolute identity "
        "difference and the rank by |log2(N50_tool/N50_real)|. R4 is "
        "excluded because its truth set is synthetic and has no real anchor.",
    )
    add_note_after_table(
        document,
        "Table 4.",
        "The real 30x ONT row is the anchor used for the delta-to-real "
        "comparison in Table 3; simulated rows are shown for context. "
        "Reference identity is matched bases divided by aligned bases from "
        "an asm5 PAF alignment of the unpolished Flye contigs against the "
        "reference. It is not a polished consensus identity and is the same "
        "metric used for every simulated row.",
    )
    add_note_after_table(
        document,
        "Table 5.",
        "This profile-matched three-seed bin ablation (E. coli and "
        "A. thaliana) is a different experiment from the six-species, "
        "two-replicate route panel reported in Table 7.",
    )
    add_note_after_table(
        document,
        "Table 6.",
        "The 5x, 10x and 15x values come from 6, 3 and 2 disjoint subsets, "
        "respectively. The 15x estimate is based on two halves and is "
        "reported as suggestive. The between-window interval is a separate, "
        "depth-independent uncertainty and is not a power statement. The "
        "subset count and per-subset coverage differ across depths, so the "
        "SD estimates are not directly comparable; the 10x estimate is "
        "based on three subsets and is not expected to be monotonic relative "
        "to the six-subset 5x estimate.",
    )
    add_note_after_table(
        document,
        "Table 7.",
        "This is the six-species B/C route panel, not the profile-matched "
        "bin ablation in Table 5. H. sapiens chr21 shows the largest "
        "directional contrast, but it is a two-replicate supporting panel, "
        "not the pre-registered two-species decision gate.",
    )


def renumber_all_text(document: Document) -> None:
    for paragraph in document.paragraphs:
        if paragraph.text.strip():
            set_paragraph_text(
                paragraph,
                renumber_figure_text(paragraph.text),
            )
    for table in document.tables:
        for row in table.rows:
            for cell in row.cells:
                for paragraph in cell.paragraphs:
                    if paragraph.text.strip():
                        set_paragraph_text(
                            paragraph,
                            renumber_figure_text(paragraph.text),
                        )


def rebuild_figure_legends(
    document: Document,
    background_captions: dict[str, str],
) -> None:
    heading = find_paragraph(document, "Figure legends")
    paragraphs = list(document.paragraphs)
    start = next(
        index for index, paragraph in enumerate(paragraphs)
        if paragraph._p is heading._p
    )
    end = len(paragraphs)
    for index in range(start + 1, len(paragraphs)):
        style_name = paragraphs[index].style.name if paragraphs[index].style else ""
        if style_name == "Heading 1":
            end = index
            break
    additions = [
        (
            "B1",
            "The panel uses the six reference genomes and the project "
            "GC-error curve. No inferential error bars are shown.",
        ),
        (
            "B2",
            "The matrix is a qualitative design summary, not a "
            "performance score. No inferential error bars are shown.",
        ),
        (
            "B3",
            "Panel a uses the frozen real E. coli ONT homopolymer table; "
            "panel b is a qualitative synthesis of the cited literature. "
            "No inferential error bars are shown.",
        ),
    ]
    anchor = heading
    for figure_number, (key, note) in enumerate(additions, start=1):
        caption = background_captions[key].strip()
        if key == "B1":
            caption = (
                "Local GC context across the six-genome panel and an "
                "illustrative real-versus-simulated GC-error curve for "
                "E. coli."
            )
            note = (
                "Panel a uses the six reference genomes; panel b uses the "
                "project GC-error curve. No inferential error bars are shown."
            )
        if not caption.endswith("."):
            caption += "."
        anchor = create_paragraph_after(
            anchor,
            f"Figure {figure_number}. {caption} {note}",
            style="Caption",
        )


def remove_page_breaks(document: Document) -> None:
    for br in list(document.element.body.iter(qn("w:br"))):
        if br.get(qn("w:type")) == "page":
            br.getparent().remove(br)


def add_line_numbers(document: Document) -> None:
    for section in document.sections:
        sect_pr = section._sectPr
        for existing in list(sect_pr.findall(qn("w:lnNumType"))):
            sect_pr.remove(existing)
        line_numbers = OxmlElement("w:lnNumType")
        line_numbers.set(qn("w:countBy"), "1")
        line_numbers.set(qn("w:restart"), "continuous")
        page_size = sect_pr.find(qn("w:pgSz"))
        if page_size is not None:
            page_size.addprevious(line_numbers)
        else:
            sect_pr.append(line_numbers)


def remove_table_shading(document: Document) -> None:
    for table in document.tables:
        for row in table.rows:
            for cell in row.cells:
                for shading in list(cell._tc.iter(qn("w:shd"))):
                    shading.getparent().remove(shading)


def set_submission_spacing(document: Document) -> None:
    for paragraph in document.paragraphs:
        style_name = paragraph.style.name if paragraph.style else ""
        if style_name in {"Title", "Heading 1", "Heading 2", "Heading 3"}:
            paragraph.paragraph_format.line_spacing = 1.0
            paragraph.paragraph_format.space_after = Pt(6)
        else:
            paragraph.paragraph_format.line_spacing = 2.0
            paragraph.paragraph_format.space_after = Pt(0)
    for table in document.tables:
        for row in table.rows:
            for cell in row.cells:
                for paragraph in cell.paragraphs:
                    paragraph.paragraph_format.line_spacing = 1.0
                    paragraph.paragraph_format.space_after = Pt(0)


def set_submission_fonts(document: Document) -> None:
    for style_name, size, bold in (
        ("Normal", 11, False),
        ("Title", 16, True),
        ("Heading 1", 14, True),
        ("Heading 2", 12, True),
        ("Heading 3", 11, True),
        ("Caption", 9, False),
    ):
        try:
            style = document.styles[style_name]
        except KeyError:
            continue
        style.font.name = "Arial"
        style.font.size = Pt(size)
        style.font.bold = bold


def convert_reference_labels(document: Document) -> None:
    paragraphs = list(document.paragraphs)
    start = next(
        index for index, paragraph in enumerate(paragraphs)
        if paragraph.text.strip() == "References"
    )
    end = len(paragraphs)
    for index in range(start + 1, len(paragraphs)):
        style_name = paragraphs[index].style.name if paragraphs[index].style else ""
        if style_name == "Heading 1":
            end = index
            break
    for paragraph in paragraphs[start + 1:end]:
        match = re.match(r"^\[(\d+)\]\s+(.*)$", paragraph.text.strip())
        if match:
            set_paragraph_text(
                paragraph,
                f"{match.group(1)}. {match.group(2)}",
            )


REFERENCE_ORDER = [
    1, 2, 3, 9, 10, 11, 12, 15, 14, 13, 16, 4, 5, 6, 7, 8,
]
REFERENCE_MAP = {
    old: new for new, old in enumerate(REFERENCE_ORDER, start=1)
}


def compress_citation_numbers(numbers: list[int]) -> str:
    ranges: list[tuple[int, int]] = []
    for number in numbers:
        if ranges and number == ranges[-1][1] + 1:
            ranges[-1] = (ranges[-1][0], number)
        elif ranges and number == ranges[-1][1]:
            continue
        else:
            ranges.append((number, number))
    parts: list[str] = []
    for start, end in ranges:
        if end - start >= 2:
            parts.append(f"{start}-{end}")
        else:
            parts.extend(str(number) for number in range(start, end + 1))
    return ", ".join(parts)


def renumber_references(document: Document) -> None:
    """Renumber citations by first mention and cite all reference entries."""
    paragraphs = list(document.paragraphs)
    reference_heading = next(
        index for index, paragraph in enumerate(paragraphs)
        if paragraph.text.strip() == "References"
    )
    old_references: dict[int, str] = {}
    reference_paragraphs = []
    for paragraph in paragraphs[reference_heading + 1:]:
        match = re.match(r"^\[(\d+)\]\s+(.*)$", paragraph.text.strip())
        if match:
            old_references[int(match.group(1))] = match.group(2)
            reference_paragraphs.append(paragraph)
    if set(old_references) != set(REFERENCE_MAP):
        raise RuntimeError("reference list is incomplete")

    citation_pattern = re.compile(r"\[(\d+(?:\s*[-,]\s*\d+)*)\]")

    def replace_group(match: re.Match[str]) -> str:
        numbers: list[int] = []
        for part in re.split(r"\s*,\s*", match.group(1)):
            if "-" in part:
                start, end = (int(value) for value in part.split("-", 1))
                numbers.extend(range(start, end + 1))
            else:
                numbers.append(int(part))
        mapped = [REFERENCE_MAP[number] for number in numbers]
        return f"[{compress_citation_numbers(mapped)}]"

    for paragraph in paragraphs[:reference_heading]:
        new_text = citation_pattern.sub(replace_group, paragraph.text)
        if new_text != paragraph.text:
            set_paragraph_text(paragraph, new_text)

    def insert_citation(prefix: str, old: str, new: str) -> None:
        for paragraph in document.paragraphs[:reference_heading]:
            if (
                paragraph.text.strip().startswith(prefix)
                and old in paragraph.text
                and new not in paragraph.text
            ):
                set_paragraph_text(
                    paragraph,
                    paragraph.text.replace(old, new, 1),
                )
                return
        raise RuntimeError(f"citation insertion failed: {old}")

    insert_citation(
        "We compared GCerrHMM with NanoSim",
        "and badread on the Level-1",
        "and badread [12] on the Level-1",
    )
    insert_citation(
        "Simulated reads were aligned",
        "same minimap2/samtools pipeline",
        "same minimap2 and samtools [13, 14] pipeline",
    )
    insert_citation(
        "Two downstream layers carry",
        "(Flye for ONT",
        "(Flye [15] for ONT",
    )
    insert_citation(
        "For phasing, real HG002",
        "processed with WhatsHap",
        "processed with WhatsHap [16]",
    )

    for index, old_number in enumerate(REFERENCE_ORDER):
        reference_text = old_references[old_number]
        if old_number == 4:
            reference_text = (
                "Badread: a read simulator for long reads. "
                "https://github.com/rrwick/Badread. "
                "Accessed 29 September 2026."
            )
        set_paragraph_text(
            reference_paragraphs[index],
            f"[{index + 1}] {reference_text}",
        )


def convert_additional_files(document: Document) -> None:
    heading = find_paragraph(document, "Supplementary material")
    set_paragraph_text(heading, "Additional files")
    replacements = {
        "Supplementary Table S1": (
            "Additional file 1. File name: Additional_file_1_species_panel"
            ".xlsx. File format: XLSX. Title: Species panel and source "
            "accessions. Description: reference accessions, genome size, "
            "GC content, GC heterogeneity and ONT run accessions."
        ),
        "Supplementary Table S2": (
            "Additional file 2. File name: Additional_file_2_level1_sub_scores"
            ".xlsx. File format: XLSX. Title: Full per-species Level-1 "
            "sub-scores and HiFi cross-platform scores. Description: "
            "replicate-level read-length, QV, GC and k-mer scores."
        ),
        "Supplementary Table S3": (
            "Additional file 3. File name: Additional_file_3_gc_bin_ablation"
            ".xlsx. File format: XLSX. Title: GC-bin ablation full table. "
            "Description: curve r and MAD by species, bin setting and seed."
        ),
        "Supplementary Table S4": (
            "Additional file 4. File name: Additional_file_4_delta_to_real"
            ".xlsx. File format: XLSX. Title: Delta-to-real per-tool table. "
            "Description: anchored-layer metrics and ranks."
        ),
        "Supplementary Table S5": (
            "Additional file 5. File name: Additional_file_5_instrument_audit"
            ".pdf. File format: PDF/DOCX. Title: Instrument audit. "
            "Description: four measurement problems, symptoms and fixes."
        ),
        "Supplementary Note 1": (
            "Additional file 6. File name: Additional_file_6_uncertainty_layers"
            ".pdf. File format: PDF/DOCX. Title: Two layers of uncertainty. "
            "Description: read-sampling and between-window components of "
            "the GC-fidelity instrument."
        ),
        "Supplementary Note 2": (
            "Additional file 7. File name: Additional_file_7_hifi_crossplatform"
            "_note.pdf. File format: PDF/DOCX. Title: HiFi cross-platform "
            "note. Description: cross-platform chain and the encoding-limited "
            "QV caveat."
        ),
    }
    for paragraph in list(document.paragraphs):
        text = paragraph.text.strip()
        for prefix, replacement in replacements.items():
            if text.startswith(prefix):
                set_paragraph_text(paragraph, replacement)
                break


def update_title_line(document: Document) -> None:
    for paragraph in document.paragraphs:
        if "Target journal: BMC Bioinformatics" in paragraph.text:
            set_paragraph_text(
                paragraph,
                "Article type: Research article | Target journal: BMC "
                "Bioinformatics | Collection \"Simulated and synthetic data "
                "in bioinformatics: statistical methods, models, and "
                "applications\" | W3.3 submission-format manuscript, "
                "2026-09-28.",
            )
            return


def update_headers(document: Document) -> None:
    for section in document.sections:
        for paragraph in section.header.paragraphs:
            if "W2 layout draft" in paragraph.text:
                set_paragraph_text(
                    paragraph,
                    "GCerrHMM - BMC Bioinformatics W3.3 submission manuscript",
                )


def copy_submission_figures(source: Path, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    for old_stem, new_stem in FIGURE_MAP.items():
        for suffix in (".png", ".pdf"):
            matches = list(source.rglob(old_stem + suffix))
            if len(matches) != 1:
                raise RuntimeError(
                    f"expected one {old_stem + suffix}, found {len(matches)}"
                )
            shutil.copy2(matches[0], output / (new_stem + suffix))


def validate_legends(document: Document) -> list[str]:
    warnings: list[str] = []
    for paragraph in document.paragraphs:
        text = paragraph.text.strip()
        if not re.match(r"^Figure \d+\.", text):
            continue
        words = len(text.split())
        if words > 300:
            warnings.append(f"{text[:40]}... has {words} words")
    return warnings


def write_extract(document: Document, output: Path) -> None:
    lines: list[str] = []
    for child in document.element.body.iterchildren():
        if child.tag.endswith("}p"):
            paragraph = Paragraph(child, document)
            text = paragraph.text.strip()
            if text:
                style_name = paragraph.style.name if paragraph.style else ""
                lines.append(f"[{style_name}] {text}" if style_name else text)
        elif child.tag.endswith("}tbl"):
            table = Table(child, document)
            for row in table.rows:
                lines.append(
                    " | ".join(cell.text.strip() for cell in row.cells)
                )
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-docx", required=True)
    parser.add_argument("--output-docx", required=True)
    parser.add_argument("--figure-source", required=True)
    parser.add_argument("--figure-output", required=True)
    args = parser.parse_args()

    source = Path(args.source_docx)
    output = Path(args.output_docx)
    document = Document(source)

    background_captions = remove_embedded_figures(document)
    missing = sorted({"B1", "B2", "B3"} - set(background_captions))
    if missing:
        raise RuntimeError(f"missing Background captions: {', '.join(missing)}")

    renumber_all_text(document)
    append_to_paragraph(
        document,
        "Widely used long-read simulators",
        "The GC context of the six-genome panel and representative "
        "GC-error curves are shown in Fig. 1.",
    )
    append_to_paragraph(
        document,
        "A simulator that conditions",
        "Fig. 2 summarizes recent simulator designs and the modelling "
        "axes they cover.",
    )
    append_to_paragraph(
        document,
        "Here we present GCerrHMM",
        "Fig. 3 places the present study within that context and "
        "evaluation chain.",
    )
    repair_incomplete_sentences(document)
    revise_scientific_text(document)
    normalize_table_labels(document)
    move_relation_section(document)
    add_table_notes(document)
    rebuild_figure_legends(document, background_captions)
    remove_page_breaks(document)
    add_line_numbers(document)
    remove_table_shading(document)
    set_submission_spacing(document)
    set_submission_fonts(document)
    renumber_references(document)
    convert_reference_labels(document)
    convert_additional_files(document)
    update_title_line(document)
    update_headers(document)

    output.parent.mkdir(parents=True, exist_ok=True)
    document.core_properties.title = (
        "GCerrHMM: a GC-aware error hidden Markov model for long-read "
        "sequencing simulation with a reproducible fidelity harness"
    )
    document.save(output)
    write_extract(document, output.with_name(output.stem + "_extract.txt"))

    figure_output = Path(args.figure_output)
    copy_submission_figures(Path(args.figure_source), figure_output)
    warnings = validate_legends(document)
    print(f"W3_SUBMISSION_DOCX={output}")
    print(f"W3_FIGURES={figure_output}")
    for warning in warnings:
        print(f"WARNING={warning}")


if __name__ == "__main__":
    main()
