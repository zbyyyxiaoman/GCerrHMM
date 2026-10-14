# Endpoint equivalence check

The second A. thaliana cohort uses the same primary endpoint as the complete
GC-error fidelity script: Pearson correlation between real and simulated
10-bin GC-stratified error curves, with 100-bp windows and a minimum-support
threshold of 20,000 aligned bases per bin.

Before scoring the second cohort, an endpoint-only implementation was run on
all six first-cohort points and compared with the complete fidelity JSON
files. Every reported value agreed exactly.

| Point | Complete fidelity | Endpoint-only | Absolute difference |
|---|---:|---:|---:|
| 1-bin, seed 44 | 0.6421581158369087 | 0.6421581158369087 | 0.0 |
| 1-bin, seed 55 | 0.6643591242646716 | 0.6643591242646716 | 0.0 |
| 1-bin, seed 66 | 0.8594777434089627 | 0.8594777434089627 | 0.0 |
| 5-bin, seed 44 | 0.9774599118175736 | 0.9774599118175736 | 0.0 |
| 5-bin, seed 55 | 0.9459997192684872 | 0.9459997192684872 | 0.0 |
| 5-bin, seed 66 | 0.9535938789478930 | 0.9535938789478930 | 0.0 |

The endpoint-only calculation omits the descriptive window-level bootstrap
and the 50-cell secondary curve. Those quantities are not part of the
pre-specified confirmation test and remain available in
`scripts/gc_error_fidelity.py`.
