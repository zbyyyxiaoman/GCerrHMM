# Final cross-tool comparison

All methods use the same reference, 10x target depth, 10,000-read training
subset, and matched-sampling Level-1 evaluation. The k-mer metric uses the
same top-N normalized sampling on both real and simulated data.

| Species | Tool | Composite | Read length | QV | GC | k-mer |
|---|---|---:|---:|---:|---:|---:|
| Ecoli | errhmm | 85.35 | 93.63 | 70.08 | 98.67 | 75.20 |
| Ecoli | NanoSim | 84.75 | 88.00 | 64.28 | 97.44 | 87.64 |
| Ecoli | PBSim3 | **93.89** | 90.99 | 100.00 | 98.66 | 80.67 |
| Ecoli | badread | 80.74 | 96.71 | 52.93 | 96.71 | 75.23 |
| Scerevisiae | errhmm | **84.77** | 95.77 | 95.87 | 82.37 | 62.98 |
| Scerevisiae | NanoSim | 76.85 | 99.71 | 87.87 | 82.56 | 29.78 |
| Scerevisiae | PBSim3 | 78.76 | 94.46 | 99.75 | 82.80 | 28.32 |
| Scerevisiae | badread | 79.93 | 97.56 | 91.50 | 83.35 | 41.16 |
| Hsapiens_chr21 | errhmm | 84.49 | 93.34 | 97.24 | 86.90 | 54.64 |
| Hsapiens_chr21 | NanoSim | 77.97 | 93.50 | 92.12 | 63.58 | 67.55 |
| Hsapiens_chr21 | PBSim3 | **84.70** | 97.25 | 100.00 | 85.16 | 50.98 |
| Hsapiens_chr21 | badread | 77.30 | 97.75 | 81.78 | 87.24 | 34.30 |

## Reading the table

- PBSim3 is the strongest distributional baseline under this composite,
  especially on QV and Ecoli.
- errHMM is strongest on Scerevisiae and competitive on Hsapiens_chr21
  once the k-mer metric is corrected.
- The result is species-dependent; it does not support a universal
  "errHMM wins" or "PBSim3 wins" claim.
- See `composite_evaluation_protocol.md` for the fairness and reporting
  caveats.
