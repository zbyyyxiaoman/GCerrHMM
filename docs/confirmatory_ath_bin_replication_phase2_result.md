# Second pre-specified A. thaliana replication cohort

Primary endpoint: paired `5-bin - 1-bin` Pearson r.

Protocol: `docs/confirmatory_ath_bin_replication_phase2_protocol.md`.

## Second cohort (seeds 77, 88, 99, 111, 122, 133)

- seed-wise deltas: `{'77': 0.3307495033952389, '88': 0.3810403144215787, '99': 0.4302173141913821, '111': 0.3387340768179975, '122': 0.2686420085565827, '133': 0.05514031096282346}`
- mean delta: `0.3007539214`
- sample SD: `0.1318378722`
- 95% t CI: `[0.1623985596, 0.4391092831]`
- paired t: `5.5878757248`, p = `0.0025322789`
- Wilcoxon p: `0.0312500000`
- pre-specified support: `True`

| Seed | 1-bin r | 5-bin r | Delta |
|---:|---:|---:|---:|
| 77 | 0.5917909261 | 0.9225404295 | 0.3307495034 |
| 88 | 0.5639706393 | 0.9450109537 | 0.3810403144 |
| 99 | 0.5329051172 | 0.9631224314 | 0.4302173142 |
| 111 | 0.6103041985 | 0.9490382753 | 0.3387340768 |
| 122 | 0.7173556069 | 0.9859976155 | 0.2686420086 |
| 133 | 0.9267526122 | 0.9818929232 | 0.0551403110 |

## First cohort (seeds 44, 55, 66)

- mean delta: `0.2370195088`
- paired t p: `0.0834162589`
- support: `False`

## Pooled 12-seed sensitivity

- mean delta: `0.2549155302`
- 95% t CI: `[0.1796121114, 0.3302189489]`
- paired t p: `0.0000127610`

## Endpoint implementation

The endpoint-only implementation uses the same reference windows, 10 equal-width GC bins, CIGAR event definition and minimum-support rule as the complete fidelity script. Before the second cohort was scored, it reproduced the six first-cohort top-level `pearson_r` values exactly; all six absolute differences were `0.0`. The unused window-level bootstrap and descriptive curve were omitted from the confirmation run and remain available in the complete fidelity script. The equivalence table is retained in `docs/confirmatory_ath_bin_replication_phase2_bundle/EQUIVALENCE.md`.

The endpoint JSON files and their SHA-256 manifest are stored under `docs/confirmatory_ath_bin_replication_phase2_bundle/results/confirm_ath_bins_phase2_20261014/status/`.

The original cross-species gate remains unchanged.
