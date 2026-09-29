# Figure captions

## Figure 1. errHMM workflow

Overview of the frozen errHMM analysis pipeline, from real ONT BAM and
reference-input preparation, through MAPQ-filtered errHMM training with
100-bp GC windows, to three simulation routes (empirical sampling,
one-bin qsHMM control, and GC-aware errHMM), followed by matched Level-1
evaluation and downstream mapping, variant-calling, and assembly analyses.

## Table 1. Species panel

Reference and read-data panel for the six analyzed genome panels. Genome
size and GC statistics are computed directly from the frozen reference
FASTA files. GC heterogeneity is the standard deviation of GC percentage
across valid non-overlapping 1-kb windows (windows with fewer than 90%
canonical A/C/G/T bases are excluded). Reference assembly accessions and
ONT read accessions identify the source data.

## Figure 2. GC-aware error-rate heatmap

Species-specific errHMM error rates by GC bin, computed from the frozen
model transition matrices. Each row is one species and each column is a
10% GC bin; values are normalized within species to their mean so that
relative GC-dependent variation is visible. This is the main GC-aware
model visualization; the row-level values are not raw absolute error
rates.

## Figure 3. Level-1 sub-score heatmap

Mean Level-1 sub-scores for six genome panels, three simulation routes,
and four evaluation dimensions. Each row represents one species-route
combination; columns show read-length similarity, QV histogram
similarity, GC distribution similarity, and 21-mer similarity. The heatmap
exposes metric-specific trade-offs that are hidden by a single composite
score. QV is profile-matched for all three routes, so the near-identical
QV column is expected to be controlled by the shared real-read quality
profile rather than by the route-specific error model.

## Figure 4. GC-aware k-mer gain

Species-wise difference in the k-mer sub-score between route C (GC-aware
errHMM) and route B (one-bin qsHMM control): route C minus route B.
Positive values indicate higher sequence-level k-mer similarity for the
GC-aware model. Bars are ordered by effect size.

## Figure 5. Cross-tool comparison

Matched-data comparison of errHMM, NanoSim, PBSim3, and badread across
E. coli, S. cerevisiae, and human chromosome 21. Facets show composite,
read-length, QV, GC, and k-mer scores. All methods use the same reference,
10x target depth, 10,000-read training subset, and matched-sampling
evaluation. PBSim3 is a sample/profile-fitted baseline; the figure should
not be interpreted as a universal biological-accuracy ranking.

## Figure 6. Coverage and state-space robustness

Left: coverage-ablation composite scores across 1x, 5x, 10x, and 30x
training coverage. Points show individual random seeds and error bars show
mean ± SD across three seeds. Right: full versus simplified state-space
models across three seeds. The figure reports the frozen robustness
analysis used for the final claims.

## Supplementary Figure S1 (working label Fig A). Aligned compound error rate

Local GC content versus aligned compound error rate with soft clips
excluded, in the frozen E. coli panel.
Each point summarizes aligned reference windows; the black line is real
ONT data and colored lines are the four compared simulators. The first
version of this panel counted soft-clipped bases as insertions and
inflated the real error rate; the corrected analysis excludes those bases
with `query_alignment_start` and `query_alignment_end`. After correction,
real and errHMM/PBSim3 error rates are of the same order (roughly 0.5-1.1%
across GC bins). The expected monotonic GC-error slope is still not
strongly reproduced, so this remains an exploratory result and should not
be used as the main evidence for GC-aware error fidelity.

## Supplementary Figure S2 (working label Fig B). Homopolymer deletion-rate curve

Deletion rate as a function of reference homopolymer length in the frozen
E. coli panel. Real ONT data show increasing deletion rate with
homopolymer length. badread shows the same qualitative direction, whereas
errHMM, NanoSim, and PBSim3 remain approximately flat. This panel does not
support a claim that errHMM reproduces homopolymer error structure; it
instead exposes a remaining simulation limitation.

## Supplementary Figure S3 (working label Fig C). k-mer composition PCA

PCA of 5-mer frequency vectors from repeated read subsamples of real ONT
reads, the three original simulation routes, and the four cross-tool
simulators. The errHMM cloud has the smallest center distance to the real
cloud in this frozen E. coli analysis, followed by route B and PBSim3.
This panel supports composition-level fidelity of errHMM but is an
exploratory two-dimensional summary rather than a formal classifier test.
badread is shown in an inset because its PC1/PC2 coordinates are on a
substantially larger scale; all eight categories use distinct tab10-style
colors.

## Supplementary Figure S4 (working label Fig D). GC heterogeneity versus k-mer gain

Species-wise scatter of reference GC heterogeneity (SD of 1-kb GC
percentage points) against the frozen route C-minus-B k-mer gain. Across
the six genome panels the linear association is weak
(Pearson r approximately -0.03). No regression line is drawn; the panel is
shown only as an exploratory Discussion scatter and is not treated as
evidence for a monotonic heterogeneity-gain relationship.

## Supplementary Figure S5 (working label Fig E). Simulator selection guide

Scenario-oriented recommendation guide rather than a performance claim.
Use PBSim3 for GC-uniform genomes when distributional similarity is the
priority; use errHMM for GC-heterogeneous eukaryotic genomes when
composition or GC fidelity is the priority; use badread when non-ideal or
dirty-read features are the target. No simulator is presented as a
universal winner.
