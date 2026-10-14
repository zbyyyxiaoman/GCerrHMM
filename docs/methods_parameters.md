# Methods parameters

## Model and training parameters

- HMM state space, full: `M`, `S`, `I`, `D1`, `D2`, `D3`, `D4+` (7 states)
- HMM state space, simplified: `M`, `S`, `I`, `D` (4 states)
- Default GC bins: 10
- Profile-matched GC-bin panels: 1, 5, 10, and 20 bins
- GC window size: 100 bp
- Minimum mapping quality for training: MAPQ >= 20
- Training-read cap in matched comparisons: 10,000 reads
- Main simulation depth: 10x unless stated otherwise
- Default random seed: 42
- GC KS sampling seed: 20260913 (local RNG; live runs are deterministic)
- Profile-matched seed list: 11, 22, 33
- Coverage seed list: 101, 202, 303
- State-space seed list: 11, 22, 33

## Level-1 composite

Sub-scores are computed on a 0-100 scale:

- read length: min(real N50, simulated N50) / max(real N50, simulated N50) * 100
- QV: (Pearson correlation of normalized QV histograms + 1) / 2 * 100
- GC: (1 - KS statistic of window GC distributions) * 100
- k-mer: matched-sampling normalized 21-mer frequency correlation * 100
- error rate (optional): 1 - absolute error-rate difference / max error rate

When available, the composite weights are:

| Component | Weight |
|---|---:|
| Read length | 0.15 |
| Error rate | 0.25 |
| QV | 0.20 |
| GC | 0.25 |
| k-mer | 0.15 |

The error-rate component was not enabled in the frozen main comparison; the
remaining weights were renormalized over the available components.

## Software versions

Versions from the Titan analysis environments:

- Minimap2: 2.31-r1302
- SAMtools: 1.23.1
- BCFtools: 1.24
- Flye: 2.9.6-b1802
- Raven assembler: 1.8.3
- NanoSim: 3.2.3
- PBSim3: 3.0.5
- badread: 0.4.2
- Python: 3.11.16
- Biopython: 1.88
- NumPy: 2.4.6
- SciPy: 1.17.1
- Matplotlib: 3.11.2

## Notes

- The k-mer comparison in the frozen route figures follows the original
  frozen Level-1 JSONs. The separate cross-tool evaluation uses the
  matched-sampling k-mer correction documented in `kmer_metric_audit.md`.
- PBSim3 is the sample/profile-fitted external baseline.
- HiFi cross-platform analysis is not reported because the available
  `Ecoli_hifi.fastq.gz` is short-read data (mean length 69.84 bp).
