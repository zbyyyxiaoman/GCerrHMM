# Manuscript source status

Current submission file:
`docs/manuscript/GCerrHMM_BMC_Research_article_W3.4_submission_20260930.docx`
with its rendered PDF and figure set in the same directory. W3.3 and earlier
files remain archived for provenance only.

The current writing skeleton is
`docs/manuscript/GCerrHMM_论文骨架_W1.3.1.docx` during the active writing
sprint. Its text extraction is `docs/manuscript/w1_3_1_extract.txt`, and its
numbers are backed by the frozen sources listed below.

The figure/table-populated layout draft is
`docs/manuscript/GCerrHMM_BMC_manuscript_W2_draft_20260927.docx` with a
rendered PDF and text extraction in the same directory. Its layout rationale
is recorded in `docs/manuscript/W2_layout_notes.md`.

The next layout revision,
`docs/manuscript/GCerrHMM_BMC_manuscript_W2_draft_20260927_v2.docx`, adds
three Background figures (B1-B3) and the refined main-figure palette and
axis styling. Its Background figure provenance is recorded in
`docs/background_figures/README.md`; the rendered PDF and text extraction
accompany the DOCX in the same directory.

The BMC submission-format conversion is
`docs/manuscript/GCerrHMM_BMC_Research_article_W3_submission_20260928.docx`.
It uses double spacing, line/page numbering, sequential figure numbering
and separate figure files. The conversion rationale is in
`docs/manuscript/W3_submission_notes.md`; the submission package is
`GCerrHMM_W3_submission_20260928/`.

The review-driven revision is
`docs/manuscript/GCerrHMM_BMC_Research_article_W3.1_submission_20260928.docx`.
It corrects the transition/emission description, one-bin interpretation,
uncertainty wording, boundary-condition framing, MAPQ rationale and table
definitions, and it names the generated Additional file 1-7 releases.
The revision package is `GCerrHMM_W3_submission_20260928_rev1/`.

The second-review revision is
`docs/manuscript/GCerrHMM_BMC_Research_article_W3.2_submission_20260929.docx`.
It shortens the abstract to 291 words, adds the real-anchor identity
definition and the depth-grid comparability note, and clarifies directional
evidence and figure-legend wording. Its response record is
`docs/manuscript/W3.2_review_response.md`.

Because the BMC submission manuscript intentionally has no embedded figures,
the figure-inclusive review copy is
`docs/manuscript/GCerrHMM_BMC_Research_article_W3.1_review_with_figures_20260928.pdf`,
and `docs/manuscript/W3.1_figures_contact_sheet.png` is a nine-figure visual
index.

The corresponding W3.2 figure-inclusive review copy is
`docs/manuscript/GCerrHMM_BMC_Research_article_W3.2_review_with_figures_20260929.pdf`.

The W3.3 style revision is
`docs/manuscript/GCerrHMM_BMC_Research_article_W3.3_submission_20260929.docx`;
the BMC expression-style audit is in
`docs/manuscript/W3_style_notes.md`.

The W3.4 numeric-provenance revision is
`docs/manuscript/GCerrHMM_BMC_Research_article_W3.4_submission_20260930.docx`.
It uses one chr21 training-alignment denominator for the insertion-rate
comparison (`3.78%`), labels the length distribution as insertion-only, and
keeps the independent remapped MAPQ-selection panel out of that sentence.

The Chinese figure-inclusive abstract is
`docs/manuscript/GCerrHMM_Chinese_figure_abstract_W3.3_20260929.docx` and
its rendered PDF is in the same directory.

The full Chinese manuscript is
`docs/manuscript/GCerrHMM_Chinese_full_manuscript_W3.3_20260929.docx`
(rendered PDF in the same directory), with its Markdown source at
`docs/manuscript/GCerrHMM_中文全文_W3.3_20260929.md`.

The reference verification and citation-order audit is
`docs/manuscript/W3_reference_audit.md`.

The `intro_draft.md`, `methods_draft.md`, and `results_r*.md` files are
historical section drafts. They remain useful for prose and method detail but
are not authoritative for numbers:

* `results_r3_alignment_draft.md` and `results_r5_assembly_phasing_draft.md`
  describe the broader main cross-tool panel and its older real-anchor text.
  They do not represent the coverage-matched chr21 30x decision panel.
* `results_r4_variant_draft.md` correctly states that the synthetic-truth
  variant layer has no real anchor and is excluded from delta-to-real.
* `results_r6_ablation_draft.md` predates the final MAD/shape separation and
  must be checked against `docs/paper_locked_conclusions.md`.

Authoritative decision sources:

```
results/framework/stats/delta_to_real_decision.md
results/framework/stats/final_layer_summary_v2.md
results/stats/gc_bins_claim_decision.md
results/stats/gc_fidelity_summary.csv
```

Authoritative 30x input contract: `docs/panel_provenance.md`.

Before quoting any historical draft, run:

```bash
PROJECT_DIR=/path/to/project bash reproduce.sh --audit-panel
```
