# Claim wording: the two branches, drafted before the last numbers land

Written 2026-09-26 04:45, while the 50-point GC curves for the four bin
settings are still being computed. Both branches are drafted now so the text
cannot drift towards whichever result arrives first.

Common ground (unaffected by either branch):

* length-aware indel encoding (single-variable ablation, E. coli +4.07 composite)
* compositional fidelity (k-mer +15.97 on E. coli v2, PCA post-analysis)
* cross-platform trainability (HiFi chain completes end to end)
* downstream delta-to-real ranking on the coverage-matched chr21 30x panel
* the instrument audit itself (MAPQ bias, insertion encoding, profile gate,
  k-mer sampling, composite limits)

---

## Branch 1 - the GC advantage survives (A. thaliana flips it)

**Headline.** GC-conditioned error modelling improves GC-conditional error
fidelity in GC-heterogeneous genomes.

**Results sentence.** On the profile-matched bin ablation, the one-bin control
had the lowest GC-fidelity correlation and the largest deviation from the real
per-GC-bin error curve in both *E. coli* and *A. thaliana* (mean r = X vs Y,
MAD = X vs Y pp; n = 3 seeds each, intervals disjoint), while the composite
summary score - which contains no GC-conditional error term - showed no
separable difference across bin settings.

**Discussion sentence.** The advantage is visible only on an instrument aligned
with the claim: the composite measures marginal distributions, so it is blind
to the difference that GC conditioning makes, and we report both outcomes with
equal prominence.

**Limitation sentence.** The cross-tool comparison does not separate the
simulators on this instrument at 10x; the effect that separates bin settings
within our model is smaller than the between-tool noise floor at that coverage.

---

## Branch 2 - contracted claim (current default)

**Headline.** A trainable, GC-aware error model with measurable GC-conditional
error structure, competitive compositional fidelity and downstream consistency.

**Abstract sentences (drop-in).**

1. We present GCerrHMM, an error HMM fitted from real ONT alignments whose
   emission model is conditioned on local GC content, and a reproducible
   harness that measures what such a model can and cannot reproduce.
2. Across six genomes the GC-aware model is directionally closer to the real
   GC-error curve than its one-bin control (+0.23 and +0.25 in correlation in
   *E. coli* and *S. cerevisiae*), but at every coverage we tested (10x, 20x,
   30x) the effect remains inside the sampling interval of the instrument, so
   we report it as directional rather than significant.
3. We quantify that noise floor: the 95 % interval on the GC-curve correlation
   narrows from W10 at 10x, to W20 at 20x, to W30 at 30x, implying that
   detecting a GC-conditional effect of this size needs approximately X depth
   or Y independent samples - a budget the field can use when designing
   simulator benchmarks.
4. The composite distributional score, which contains no GC-conditional error
   term, does not separate the two models; we therefore report it as a
   descriptive summary and name the instrument each claim rests on.

**Results sentences.**

* *Instrument alignment.* "The composite aggregates read-length, QV, GC-content
  and k-mer terms; none of them is a function of local GC content, so it cannot
  test the property the method changes. We therefore evaluate GC conditioning
  on a GC-stratified error curve with an explicit interval."
* *Directional consistency.* "The GC-aware model leads its one-bin control on
  the curve correlation in both replicates of both species tested (E. coli
  +0.176 vs -0.054; S. cerevisiae +0.144 vs -0.101 at 20x), and in the
  cross-tool panel no simulator - including GCerrHMM - separates from zero at
  10x."
* *Noise floor.* "Interval width as a function of coverage: 10x W10, 20x W20,
  30x W30 (Figure X). The effect size we observe is below the 10x floor, which
  is why the replicated comparison is reported as not discriminable rather
  than as a win."
* *Independent assets.* indel encoding ablation, compositional fidelity,
  HiFi trainability, downstream delta-to-real ranking.

**Limitation sentences (verbatim candidates).**

* The Level-1 composite is a distributional similarity score and contains no
  GC-conditional error term; conclusions about GC conditioning rest on the
  stratified-error instrument and the controlled bin ablation.
* The k-mer sub-score correlates frequency vectors over the intersection of the
  per-side top-k maps; for genomes with a small intersection
  (*D. melanogaster*: r = 0.12-0.35) it has little discriminating power, and
  its replicate-to-replicate spread within one route reaches 16.7 points on a
  0-100 scale.
* With 1000 reads per side and a paired bootstrap, 11 of 12 route/replicate
  comparisons of the six-species panel have a difference interval containing
  zero, so that panel is reported as not discriminable.
* At 10-30x, the GC-conditional effect we can generate is smaller than the
  noise floor of the fidelity instrument; we quantify the required depth
  rather than claim the effect.

---

## Shared "instrument audit" paragraph (both branches)

Four measurement problems were found and fixed during development, each of
which had produced plausible but wrong numbers: a MAPQ threshold that removed
the reads carrying most of the error signal; a single-insertion-state encoding
that capped every insertion at 1 bp; a missing read-length profile that
silently fell back to a platform default; and a k-mer sub-score whose
replicate-to-replicate spread exceeded the between-route differences it was
used to rank. We report them because they determine how any simulator
benchmark in this space should be read, and because the pipeline now fails
loudly instead of degrading quietly.
