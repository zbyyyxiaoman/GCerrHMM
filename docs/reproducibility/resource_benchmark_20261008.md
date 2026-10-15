# GCerrHMM resource benchmark (2026-10-08)

Measurement host: Titan (Ubuntu 24.04.1 LTS, Linux kernel
6.8.0-146-generic; AMD Ryzen 9 3950X 16-core/32-thread CPU; 125 GiB RAM),
`nice -n 19`, GNU `/usr/bin/time -v`.

Inputs: E. coli reference, real ONT aligned BAM, and the matching real-read
profile. Training used 10 GC bins, 100-bp windows, MAPQ >= 20 and at most
20,000 reads. Simulation used 10x coverage and seed 42.

| Stage | Workers | Wall seconds | Peak RSS MB | Reads | Reads/s | Speedup |
|---|---:|---:|---:|---:|---:|---:|
| Train GC-HMM | 1 | 88.12 | 72.73 | - | - | 1.00x |
| Simulate 10x | 1 | 29.56 | 106.73 | 4665 | 157.8 | 1.00x |
| Simulate 10x | 2 | 17.11 | 122.40 | 4665 | 272.6 | 1.73x |
| Simulate 10x | 4 | 10.64 | 123.02 | 4665 | 438.4 | 2.78x |

The one-, two- and four-worker FASTQ files were SHA-256-identical:

`f31abd7c5607ac73b4a9fade0c837850b63d13136065f8d6b8fa33577c856b69`

The associated machine-readable table is
`docs/reproducibility/resource_benchmark_20261008.csv`.
