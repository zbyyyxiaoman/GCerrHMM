# GC-aware improvement reproduction

- Species: `Ecoli`
- Coverage: `10.0x`
- Seed: `42`
- Simulation reference: `data/truth/Ecoli_variant_ref.fa`

| Route | Composite | Read length | QV | GC | k-mer |
|---|---:|---:|---:|---:|---:|
| sample_A | 85.31 | 99.24 | 70.09 | 97.90 | 70.72 |
| errhmm_1bin | 86.69 | 99.19 | 70.10 | 98.68 | 76.32 |
| errhmm_gc | 86.92 | 99.19 | 70.10 | 98.73 | 77.40 |

## GC-aware minus one-bin control

| Metric | Delta |
|---|---:|
| composite | +0.23 |
| read_length | +0.00 |
| qv | +0.00 |
| gc | +0.05 |
| kmer | +1.08 |

The gate is reported, not forced. A failed gate is preserved as a
negative result and must not be rewritten into a positive claim.
