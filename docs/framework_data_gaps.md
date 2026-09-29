# Framework data gaps

Date: 2026-09-22

## Registered real HiFi sources

The source registry is `config/hifi_sources.json`. Each row is an ENA/SRA
run with an accession, official file size, official MD5, instrument model,
and library name. The validator writes both JSON and Markdown QC output.

| Species | Accession | Instrument | Source FASTQ size |
|---|---|---:|---:|
| E. coli | ERR12723508 | Sequel IIe | 0.88 GB |
| S. cerevisiae | SRR31637145 | Sequel II | 6.90 GB |
| A. thaliana | SRR14728885 | Sequel II | 18.53 GB |
| D. melanogaster | SRR9969842 | Sequel II | 4.70 GB |
| M. musculus | SRR28703375 | Sequel II | 40.18 GB |
| H. sapiens | ERR13110527 | Sequel II | 32.07 GB |

The downloader uses `aria2c` with resume support, verifies exact byte size
and MD5, and records the local manifest. The validator then measures read
count, total bases, mean/median length, and mean QV from the actual FASTQ.

## Remaining framework gaps

1. **External downstream results.** The frozen paper tables currently only
   contain read-level external comparisons. The framework runner fills
   mapping, variant, and assembly cells for E. coli, S. cerevisiae, and
   human chr21.
   The first external rerun exposed a fairness bug: those reads had been
   generated from the original reference, while the internal routes were
   generated from the truth-variant reference. The corrected branch
   rebuilds all external tools from `*_variant_ref.fa` and writes `_fair`
   outputs; the original files remain untouched as provenance.
2. **Native PBSim3-HMM.** The frozen `PBSim3` rows are sample-mode results.
   They must not be described as the submitted errHMM-GC model or as the
   PBSim3 internal HMM mode.
3. **Phasing.** The current simulated panels are haploid-looking in the
   analysis sense: all reads are generated from one simulated donor
   sequence. A switch-error claim requires a diploid source, a phased truth
   VCF, and a phasing-aware evaluation. This is intentionally not inferred
   from the existing panels.
4. **Structural variants.** The current variant metric is SNP/indel F1
   from `bcftools call`. An SV layer requires a separate truth set and an
   SV caller; it is not fabricated from the existing VCF metrics.

## Execution entry points

```bash
# Official HiFi data (inside Titan tmux)
SPECIES=Ecoli,Scerevisiae,Hsapiens_chr21 \
  bash scripts/titan_fetch_hifi.sh

# Four-layer external evaluation (inside Titan tmux)
SPECIES=Ecoli,Scerevisiae,Hsapiens_chr21 \
TOOLS=errhmm,nanosim,pbsim,badread \
PARALLEL=4 THREADS=4 \
  bash scripts/titan_framework_external.sh

# Status
bash scripts/titan_framework_status.sh
```

All compute commands use `nice -n 19`, four concurrent tasks, and four
threads per task, for a total concurrency ceiling of 16 threads.
