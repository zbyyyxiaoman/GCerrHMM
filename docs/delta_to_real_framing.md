# Evaluation framing: delta-to-real

Adopted 2026-09-23, **before** the corrected numbers exist, so that the
criterion is fixed ahead of the results rather than fitted to them.

## Why the framing changed

Every published long-read simulator here - including GCerrHMM - is trained
from aligned reads. That means all of them inherit a structural blind spot:
the training signal only contains bases that survived alignment.

Measured on the human chr21 real read set (8,189 reads, 6.62x):

| | real ONT | GCerrHMM simulation |
|---|---|---|
| base identity (same mapper, same code) | 0.8663 | 0.9679 |
| bases inside the primary alignment | ~50 % | ~96 % |
| reads with MAPQ < 20 | 58.0 % | - |

The simulator is therefore *cleaner than the data it is meant to imitate*.
Under an "absolute score, higher is better" reading, GCerrHMM wins the
alignment and variant layers. Under a "closer to real is better" reading, the
ordering can reverse: NanoSim's chr21 identity of 0.9020 is 0.036 away from
the real value of 0.8663, while GCerrHMM's 0.9679 is 0.102 away.

Both readings are defensible; only one can be the headline. We adopt
**delta-to-real**.

## Consequences

1. **All "leading" language is withdrawn** until re-measured under this
   framing. The cross-tool ranking itself stays reportable, because every
   simulator was created and evaluated under an identical protocol; what is
   withdrawn is the claim that GCerrHMM is *better* for being higher on
   absolute read-fidelity metrics.
2. **Every applicable layer carries a real-data anchor measured with the same
   pipeline**: real ONT reads through the same mapper (R3), the same caller
   (R4) and the same assembler (R5). Where the real data cannot support a
   layer - the 6.62x ONT set cannot be assembled by Flye - the anchor is
   reported as unavailable rather than substituted with a different platform.
3. **A model that scores closer to the real anchor is not automatically
   better either.** A simulator can land near the real value for the wrong
   reason (for example by producing uniformly noisy reads rather than the real
   mixture of clean and bad reads). The GC-stratified error curve and the
   ability to condition on GC therefore remain the discriminating evidence,
   and they are re-drawn under this framing (Fig A/B re-analysis).
4. **The structural blind spot is a stated limitation, not a defect of this
   work.** Alignment-based training cannot see unmapped reads, soft-clipped
   tails, or low-mappability regions. NanoSim and PBSim3-errhmm share it.
   Discussion states the boundary explicitly.

## What replaces the old headline

Old: "GCerrHMM leads on alignment, variant calling and assembly for human
chr21."

New, pending re-measurement: "GCerrHMM reproduces the GC-dependent structure
of the error rate that sampling-based simulators miss; under a delta-to-real
reading its absolute accuracy on human chr21 is X away from the real anchor,
against Y for the strongest baseline. Its residual gap to real data comes from
the alignment-based training blind spot shared by all simulators in this
panel."

The replacement claim is weaker but survives review, and the GC-conditioning
result is unaffected by the framing change because it is a controlled
comparison between two models trained from the same reads.

## Pre-registered decision rule (recorded 2026-09-25, before the 30x anchors exist)

The rule is written down **before** the coverage-matched real anchor has been
read, so that the headline cannot be chosen after seeing which way the numbers
fall.

Operational definitions, all on the coverage-matched chr21 30x panel (four
tools: GCerrHMM, NanoSim, PBSim3-sample, badread; one real anchor measured by
the identical pipeline):

| layer | quantity | distance to real |
|---|---|---|
| R3 alignment | `base_identity` | \|tool - real\| |
| R4 variant | `f1_score` (SNP and indel reported alongside) | \|tool - real\| |
| R5 assembly | `reference_identity` and `n50` | Borda rank over the two sub-metrics (rank by \|Δidentity\| plus rank by \|log2(N50_tool / N50_real)\|, averaged) |

A tool "places top-2" in a layer when its distance ranks first or second among
the four tools.

**Rule.** If GCerrHMM places top-2 in at least two of the three layers, the
headline is *downstream consistency on a GC-heterogeneous eukaryotic genome*.
Otherwise the headline falls back to *compositional fidelity on GC-heterogeneous
genomes plus a trainable framework*, and the Discussion states explicitly that
PBSim3-sample's direct sampling reproduces the real error spectrum more
faithfully at 30x.

Two further commitments:

1. The layer in which GCerrHMM does **not** place top-2 is reported in the
   abstract-level narrative, not buried in the supplement.
2. If the two branches disagree between their primary metric and their
   sub-metrics (for example a top-2 overall F1 driven by SNPs while the indel
   F1 is last), that split is reported as the result rather than resolved in
   favour of the friendlier reading.
