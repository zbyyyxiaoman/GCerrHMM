# R1 provenance

Every score is recomputed from the frozen Level-1 sub-score fields in the listed JSON file. The normalized locked weights are `0.15 read_length + 0.20 qv + 0.25 gc + 0.15 kmer`.

| Species | Route | Replicate | Frozen file | Composite line | Sub-scores line | File composite | Locked recompute | Delta |
|---|---|---:|---|---:|---:|---:|---:|---:|
| E. coli | A: Sample | 1 | `results\tables\level1_Ecoli_route_A_sample_r1.json` | 2711 | 2713 | 80.370000 | 80.367333 | 0.002667 |
| E. coli | A: Sample | 2 | `results\tables\level1_Ecoli_route_A_sample_r2.json` | 2711 | 2713 | 81.440000 | 81.443333 | 0.003333 |
| E. coli | B: qsHMM (1-bin) | 1 | `results\tables\level1_Ecoli_route_B_qshmm_r1.json` | 2711 | 2713 | 85.260000 | 85.258667 | 0.001333 |
| E. coli | B: qsHMM (1-bin) | 2 | `results\tables\level1_Ecoli_route_B_qshmm_r2.json` | 2711 | 2713 | 88.940000 | 88.938667 | 0.001333 |
| E. coli | C: errHMM (GC-aware) | 1 | `results\tables\level1_Ecoli_route_C_errhmm_r1.json` | 2711 | 2713 | 86.920000 | 86.921333 | 0.001333 |
| E. coli | C: errHMM (GC-aware) | 2 | `results\tables\level1_Ecoli_route_C_errhmm_r2.json` | 2711 | 2713 | 87.140000 | 87.135333 | 0.004667 |
| S. cerevisiae | A: Sample | 1 | `results\tables\level1_Scerevisiae_route_A_sample_r1.json` | 2711 | 2713 | 84.440000 | 84.444667 | 0.004667 |
| S. cerevisiae | A: Sample | 2 | `results\tables\level1_Scerevisiae_route_A_sample_r2.json` | 2711 | 2713 | 76.530000 | 76.530667 | 0.000667 |
| S. cerevisiae | B: qsHMM (1-bin) | 1 | `results\tables\level1_Scerevisiae_route_B_qshmm_r1.json` | 2711 | 2713 | 86.510000 | 86.510000 | 0.000000 |
| S. cerevisiae | B: qsHMM (1-bin) | 2 | `results\tables\level1_Scerevisiae_route_B_qshmm_r2.json` | 2711 | 2713 | 79.120000 | 79.119333 | 0.000667 |
| S. cerevisiae | C: errHMM (GC-aware) | 1 | `results\tables\level1_Scerevisiae_route_C_errhmm_r1.json` | 2711 | 2713 | 81.310000 | 81.314000 | 0.004000 |
| S. cerevisiae | C: errHMM (GC-aware) | 2 | `results\tables\level1_Scerevisiae_route_C_errhmm_r2.json` | 2711 | 2713 | 89.620000 | 89.616667 | 0.003333 |
| A. thaliana | A: Sample | 1 | `results\tables\level1_Athaliana_route_A_sample_r1.json` | 2711 | 2713 | 97.430000 | 97.435333 | 0.005333 |
| A. thaliana | A: Sample | 2 | `results\tables\level1_Athaliana_route_A_sample_r2.json` | 2711 | 2713 | 95.200000 | 95.197333 | 0.002667 |
| A. thaliana | B: qsHMM (1-bin) | 1 | `results\tables\level1_Athaliana_route_B_qshmm_r1.json` | 2711 | 2713 | 93.310000 | 93.316000 | 0.006000 |
| A. thaliana | B: qsHMM (1-bin) | 2 | `results\tables\level1_Athaliana_route_B_qshmm_r2.json` | 2711 | 2713 | 96.360000 | 96.364667 | 0.004667 |
| A. thaliana | C: errHMM (GC-aware) | 1 | `results\tables\level1_Athaliana_route_C_errhmm_r1.json` | 2711 | 2713 | 95.070000 | 95.073333 | 0.003333 |
| A. thaliana | C: errHMM (GC-aware) | 2 | `results\tables\level1_Athaliana_route_C_errhmm_r2.json` | 2711 | 2713 | 93.800000 | 93.794000 | 0.006000 |
| D. melanogaster | A: Sample | 1 | `results\tables\level1_Dmelanogaster_route_A_sample_r1.json` | 2711 | 2713 | 85.870000 | 85.872667 | 0.002667 |
| D. melanogaster | A: Sample | 2 | `results\tables\level1_Dmelanogaster_route_A_sample_r2.json` | 2711 | 2713 | 73.950000 | 73.952667 | 0.002667 |
| D. melanogaster | B: qsHMM (1-bin) | 1 | `results\tables\level1_Dmelanogaster_route_B_qshmm_r1.json` | 2711 | 2713 | 85.120000 | 85.120000 | 0.000000 |
| D. melanogaster | B: qsHMM (1-bin) | 2 | `results\tables\level1_Dmelanogaster_route_B_qshmm_r2.json` | 2711 | 2713 | 84.090000 | 84.094000 | 0.004000 |
| D. melanogaster | C: errHMM (GC-aware) | 1 | `results\tables\level1_Dmelanogaster_route_C_errhmm_r1.json` | 2711 | 2713 | 82.720000 | 82.718000 | 0.002000 |
| D. melanogaster | C: errHMM (GC-aware) | 2 | `results\tables\level1_Dmelanogaster_route_C_errhmm_r2.json` | 2711 | 2713 | 92.530000 | 92.527333 | 0.002667 |
| M. musculus chr19 | A: Sample | 1 | `results\tables\level1_Mmusculus_chr19_route_A_sample_r1.json` | 2711 | 2713 | 91.420000 | 91.421333 | 0.001333 |
| M. musculus chr19 | A: Sample | 2 | `results\tables\level1_Mmusculus_chr19_route_A_sample_r2.json` | 2711 | 2713 | 92.630000 | 92.627333 | 0.002667 |
| M. musculus chr19 | B: qsHMM (1-bin) | 1 | `results\tables\level1_Mmusculus_chr19_route_B_qshmm_r1.json` | 2711 | 2713 | 91.720000 | 91.724667 | 0.004667 |
| M. musculus chr19 | B: qsHMM (1-bin) | 2 | `results\tables\level1_Mmusculus_chr19_route_B_qshmm_r2.json` | 2711 | 2713 | 91.100000 | 91.095333 | 0.004667 |
| M. musculus chr19 | C: errHMM (GC-aware) | 1 | `results\tables\level1_Mmusculus_chr19_route_C_errhmm_r1.json` | 2711 | 2713 | 91.520000 | 91.525333 | 0.005333 |
| M. musculus chr19 | C: errHMM (GC-aware) | 2 | `results\tables\level1_Mmusculus_chr19_route_C_errhmm_r2.json` | 2711 | 2713 | 91.630000 | 91.627333 | 0.002667 |
| H. sapiens chr21 | A: Sample | 1 | `results\tables\level1_Hsapiens_chr21_route_A_sample_r1.json` | 2711 | 2713 | 83.560000 | 83.560667 | 0.000667 |
| H. sapiens chr21 | A: Sample | 2 | `results\tables\level1_Hsapiens_chr21_route_A_sample_r2.json` | 2711 | 2713 | 81.340000 | 81.335333 | 0.004667 |
| H. sapiens chr21 | B: qsHMM (1-bin) | 1 | `results\tables\level1_Hsapiens_chr21_route_B_qshmm_r1.json` | 2711 | 2713 | 78.680000 | 78.685333 | 0.005333 |
| H. sapiens chr21 | B: qsHMM (1-bin) | 2 | `results\tables\level1_Hsapiens_chr21_route_B_qshmm_r2.json` | 2711 | 2713 | 83.500000 | 83.505333 | 0.005333 |
| H. sapiens chr21 | C: errHMM (GC-aware) | 1 | `results\tables\level1_Hsapiens_chr21_route_C_errhmm_r1.json` | 2711 | 2713 | 81.210000 | 81.213333 | 0.003333 |
| H. sapiens chr21 | C: errHMM (GC-aware) | 2 | `results\tables\level1_Hsapiens_chr21_route_C_errhmm_r2.json` | 2711 | 2713 | 83.510000 | 83.516000 | 0.006000 |

Maximum file-vs-recompute delta: `0.00600000`. The same sub-score fields are used by fig3; the k-mer C-minus-B values used by fig4 are recomputed from the same frozen `kmer` sub-scores.
