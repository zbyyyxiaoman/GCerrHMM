# GCerrHMM main-figure production status

Date: 2026-09-22 (historical)

> **Final status (2026-09-27 17:10):** the table below records the earlier
> production backlog and is retained as provenance. The final six-figure
> snapshot is `docs/gcerrhmm_main_figures_final_20260927_1704/`; its plotted
> heatmaps pass a strict finite-cell gate, and all hashes are appended to
> `docs/data_freeze_manifest.md`.

This file follows the teacher's six-section main-figure plan. The
production generator is `scripts/build_gcerrhmm_main_figures.py`.

| Main figure | Panel content | Available now | Missing before final submission |
|---|---|---|---|
| Figure 1 | dual-platform workflow, GC-state mechanism, species GC panel | ONT/HiFi workflow, trained GC matrices, species table | final platform wording and styling |
| Figure 2 | read length, QV, GC, k-mer; tool x species; ONT/HiFi facets | ONT matched panel | HiFi matched simulation and external HiFi baselines |
| Figure 3 | mapping rate, base identity, reads identity, MAPQ | mapping rate and MAPQ | base-level identity pass and HiFi alignment panel |
| Figure 4 | SNP, indel, SV precision/recall/F1 | combined SNP/indel F1, precision, recall on ONT | split SNP/indel metrics, SV caller and truth, HiFi panel |
| Figure 5 | assembly N50/NG50/identity/BUSCO; phase block N50/switch error | legacy ONT assembly metrics | hifiasm/Flye rerun, BUSCO, whatshap diploid chr21 |
| Figure 6 | GC-bin, coverage, state-space ablations | complete three-seed frozen panel | final GCerrHMM naming and style |
| Supplementary scaling | process-level parallel simulation | E. coli 10x: 4 workers = 3.22x speedup, identical read/base statistics | larger-genome scaling run |

## Platform policy

- HiFi assembly uses `hifiasm`.
- ONT assembly uses `Flye`.
- The current legacy assembly values were generated with the historical
  raven path and are retained only until the primary-assembler rerun.
- Phasing uses `whatshap`; the truth source is GIAB HG002 benchmark VCF.
  `config/phasing_sources.json` records the HG002 HiFi run and the
  GRCh38 benchmark VCF.

## Reproduction commands

```bash
# core algorithmic improvement
SPECIES=Ecoli COVERAGE=10 SEED=42 THREADS=4 bash reproduce.sh --gc-demo

# framework tables and all main figures
bash reproduce.sh --framework-figures

# real HG002 truth (small file; run in tmux on Titan)
FETCH_HIFI=0 bash scripts/titan_fetch_phasing.sh

# optional full HG002 HiFi download for the phasing panel
FETCH_HIFI=1 bash scripts/titan_fetch_phasing.sh
```

No panel is converted from a missing measurement into a zero value.

## QV audit note

The E. coli HiFi FASTQ (`ERR12723508`) was checked for the suspected
off-by-33 error. Raw ASCII quality mean is `113.28`; Phred after
subtracting 33 is `80.28`; `pysam.get_quality_array()` reports `80.08`.
The absolute QV `83.45` is therefore not a script encoding bug, but the
source FASTQ uses many `~` (Q93) quality characters. Its QV distribution
is treated as encoding-limited for formal HiFi claims until a
quality-complete BAM/source is used.

## Human HiFi critical path

The ENA BAM for the current human HiFi run is unaligned (`FLAG=4`, no
`@SQ` records), so indexed remote extraction of chr21 from that BAM is
not possible. The working alternative is the public GIAB S3 aligned BAM
`HG002.Sequel.15kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam`;
`samtools view <S3-BAM> 21 | samtools fastq` streams only chr21 reads
without downloading the 67 GB BAM. The extraction is running in
`hg002_chr21` and the post-extract watcher `hg002_post` will align the
reads to GRCh38 chr21 and run whatshap/hifiasm. No full-genome result is
claimed before that pipeline completes.
