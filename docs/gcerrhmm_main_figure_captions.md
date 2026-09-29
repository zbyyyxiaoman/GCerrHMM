# GCerrHMM main-figure captions

## Figure 1. GCerrHMM overview and innovation

(a) Dual-platform workflow from aligned ONT or PacBio HiFi BAM through
CIGAR-derived error-state learning, GC-conditioned transition modeling,
matched simulation, and four-layer evaluation. (b) Model-implied error
state varies across GC bins for representative genomes. (c) Species panel
with the analyzed reference accession, genome size, reference GC content,
1-kb GC heterogeneity, and real ONT run accession.

## Figure 2. Reads-level cross-tool comparison

Matched-data comparison of GCerrHMM, NanoSim, PBSim3-sample, PBSim3-errhmm,
and badread using read-length, QV, GC, and k-mer similarity. The current
measured facet is ONT; the HiFi facet is emitted only after the same
training, simulation, and matched-sampling evaluation has completed. All
plotted read-level cells are measured; PBSim3-errhmm is retained because
these read statistics do not depend on its Q0 quality strings.

## Figure 3. Alignment-level comparison

Mapping rate, mean MAPQ, and base-level identity after a common
minimap2/samtools alignment pipeline. All plotted cells are measured.
PBSim3-errhmm is retained here because its mapping statistics are valid,
even though its Q0 placeholder qualities make standard variant calling
non-evaluable. The H. sapiens chr21 row uses the coverage-matched 30x panel.

## Figure 4. Variant-calling comparison

Combined SNP/indel F1 for E. coli, S. cerevisiae, and H. sapiens chr21
against the same truth VCF, plus SNP-only, indel-only, and SV spike-in F1
for the two species where every plotted tool has that measurement.
The H. sapiens chr21 row uses the coverage-matched 30x panel; this layer is
reported as a cross-tool comparison and is excluded from delta-to-real.
PBSim3-errhmm is omitted because its Q0 placeholder qualities make the
standard caller non-evaluable. H. sapiens chr21 is shown in the combined-F1
panel only because PBSim3-sample has no SNP/indel split.

## Figure 5. Assembly and phasing

Assembly N50 and reference identity for the measured panel, together with
the phasing status. HiFi is assembled with hifiasm, ONT with Flye, and
diploid phasing is evaluated with whatshap against a phased truth VCF.
Historical raven values are not relabeled as hifiasm or Flye results, and
PBSim3-errhmm is omitted because no assembly run exists for it.
The H. sapiens chr21 row uses the coverage-matched 30x panel; N50 colors
span 0-40 Mb so that the real-anchor and simulated ranges remain visible.

## Figure 6. Ablation and robustness

(a) Profile-matched GC-bin ablation for E. coli and A. thaliana across
three seeds. (b) Training-coverage ablation from 1x to 30x. (c) Full
versus simplified state-space models. Error bars are mean ± SD.

## Supplementary parallel-scaling panel

Process-level parallel read simulation measured with one, two, and four
workers at 10x E. coli coverage. The figure reports wall time and speedup
against the one-worker baseline. Read count, total bases, mean length, and
N50 are checked to remain identical across worker counts.
