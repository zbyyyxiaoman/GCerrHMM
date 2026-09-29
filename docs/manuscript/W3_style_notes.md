# BMC Bioinformatics style review notes

This note records style patterns learned from published BMC Bioinformatics
method, benchmark and simulator papers. The source papers were read only for
structure and rhetorical strategy; no sentence was copied into the
manuscript.

## Source papers

| Paper | DOI | Useful style lesson |
|---|---|---|
| MitoHiFi, BMC Bioinformatics 2023 | 10.1186/s12859-023-05385-y | A tool paper opens Background with the enabling technology and the missing capability, then uses Results subsection headings that name the actual operation being evaluated. |
| Long-read shotgun metagenomics benchmark, BMC Bioinformatics 2022 | 10.1186/s12859-022-05103-0 | A benchmark paper states the evaluation gap early, then organises Results by endpoint instead of by tool. |
| Statistical power for cluster analysis, BMC Bioinformatics 2022 | 10.1186/s12859-022-04675-1 | Statistical papers separate simulation design, assumptions, effect size and uncertainty, then state what each uncertainty component can and cannot answer. |
| Decided sample size in machine-learning applications, BMC Bioinformatics 2023 | 10.1186/s12859-023-05156-9 | Sample-size and effect-size limitations are stated explicitly and early, not hidden in a final caveat. |
| MOV&RSim, BMC Bioinformatics 2025 | 10.1186/s12859-025-06292-0 | A recent simulator paper states the contribution in a compact "this work" paragraph, compares against existing tools, and closes with a distinct limitations subsection. |
| vcfsim, BMC Bioinformatics 2026 | 10.1186/s12859-026-06453-9 | A short simulation-tool paper uses a direct abstract, Implementation/Validation/Limitations headings, and an Availability and requirements block. |

## Patterns applied to W3.3

* **Direct gap-first Background.** The first Background paragraph now moves
  from the benchmarking requirement to the missing measurement problem, then
  introduces the tool.
* **Outcome-first Results.** Results paragraphs state the comparison and the
  measured direction first; methodological justification follows.
* **Uncertainty in plain language.** Read-sampling uncertainty and
  between-window uncertainty are not treated as interchangeable. The paper
  explicitly states which component is a power statement and which is not.
* **Directional rather than rhetorical significance.** The paper reports
  effect direction, uncertainty width and the pre-registered gate result in
  the same paragraph rather than using vague phrases such as "trend toward
  significance".
* **Human case highlighted without promotion.** The strongest directional
  signal, H. sapiens chr21, is named explicitly, followed immediately by its
  status as a two-replicate supporting panel rather than part of the
  pre-registered decision gate.
* **Short, scannable Discussion opening.** The Discussion begins with the
  question and the answer in four sentences, then moves to interpretation,
  limitations and relation to prior tools.
* **Section-level rather than claim-level repetition.** The composite,
  stratified instrument, downstream consistency and pre-registered decision
  are each discussed once in their own subsection and cross-referenced
  elsewhere.

## Deliberately avoided

* No exact sentences or distinctive phrases from the source papers were
  reused.
* No promotional claims of universal superiority were added.
* No new statistical test was introduced merely to make the human result
  sound significant.
* No author, funding, repository or DOI facts were fabricated.
