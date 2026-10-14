# Pre-specified second A. thaliana replication cohort

Status: frozen before the new cohort is launched.

## Purpose

The first independent A. thaliana cohort (seeds 44, 55 and 66) was
directionally positive but did not meet the pre-specified p < 0.05
threshold. This protocol defines a second, larger independent cohort for the
same single primary contrast. It is a new experiment, not a retrospective
change to the original cross-species gate.

## Fixed design

* Species: *Arabidopsis thaliana*.
* Reference, real-read alignment, read-length/QV profile and trained models:
  the frozen v2 inputs already used by the panel.
* Coverage: 10x, profile-matched.
* Bin counts: 1 and 5.
* New seeds: 77, 88, 99, 111, 122 and 133.
* Pairing: each seed is simulated once under each bin count; the two
  GC-curve correlations are compared within seed.
* Resource ceiling: four workers per point, four points in parallel
  (16 process workers total), `nice -n 19`, tmux.
* Run tag: `phase2_20261014`.

## Primary endpoint

`delta_r = Pearson r(5-bin GC-error curve) - Pearson r(1-bin GC-error curve)`

for the six new seeds. The primary test is a two-sided paired t-test at
alpha = 0.05. The result is reported with the mean delta, seed-wise deltas,
sample SD and a 95% t-based confidence interval.

The cohort is considered confirmatory only if:

1. mean `delta_r > 0`; and
2. the two-sided paired t-test p-value is below 0.05.

Both conditions are fixed before the run. No seed, endpoint or bin count can
be changed after seeing the result.

## Secondary analysis

After the primary analysis, a pooled paired analysis over all 12 seeds
(11, 22, 33, 44, 55, 66, 77, 88, 99, 111, 122, 133) is reported as a
sensitivity analysis. The original six-species gate remains unchanged and
is not reopened by this cohort.

## Outputs

Source files are written to
`results/confirm_ath_phase2_20261014/status/` on Titan. The chain writes
`CONFIRM_ATH_ALL_POINTS_DONE` only after all 12 points finish. No GitHub or
Zenodo upload is part of this protocol.
