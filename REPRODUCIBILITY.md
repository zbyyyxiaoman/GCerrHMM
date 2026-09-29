# Reproducibility and data availability

## Environment

```bash
conda env create -f environment.yml
conda activate errhmm
```

The external tool comparison environment is separate:

```bash
conda env create -f environment-tools.yml
conda activate toolcompare
```

Python package versions are also listed in `requirements.txt`. The complete
software version record is in `docs/methods_parameters.md`.

## Frozen artifact smoke check

The repository tracks the final publication artifacts and checksum files.
Run:

```bash
bash reproduce.sh --smoke
```

This verifies:

- the R1 Level-1 overall CSV/Markdown/provenance files;
- the frozen main and innovation figure archives;
- the recorded SHA-256 values without rerunning simulations.

To rebuild the R1 table from frozen Level-1 JSON files when the full local
result tree is available:

```bash
bash reproduce.sh --r1
```

To verify every checksum recorded in `docs/data_freeze_manifest.md`:

```bash
PROJECT_DIR=/path/to/project bash reproduce.sh --full-check
```

The manifest is append-only. Verification uses the latest recorded hash for
each relative artifact and ignores superseded entries in older freeze blocks.
This mode requires the downloaded data and frozen result tree; a bare code
checkout should use `--smoke`, `--tests` and `--static-check` instead.

## Core algorithmic improvement reproduction

The central algorithmic claim is not a static figure. It is the
GC-conditioned transition model learned from a real aligned BAM. The
one-command reproduction loop trains both the GC-aware model and the
one-bin control from the same BAM, simulates matched reads with identical
coverage and seed, evaluates both routes plus the empirical-sampling
baseline, and writes the C-minus-B delta:

```bash
SPECIES=Ecoli COVERAGE=10 SEED=42 THREADS=16 JOBS=3 bash reproduce.sh --gc-demo
```

For a fast deterministic smoke run of the same chain:

```bash
SPECIES=Ecoli SEED=42 THREADS=16 JOBS=3 bash reproduce.sh --quick-reproduce
```

`--quick-reproduce` downloads the public E. coli reference and ONT FASTQ,
checks the configured FASTQ checksum, builds the aligned BAM with
minimap2/samtools, and then runs the same deterministic quick chain. If the
inputs already exist, use `--quick-demo` directly.

The command refuses to overwrite an existing run directory. Outputs are
written under `results/gc_improvement/<species>_<timestamp>/` and include:

- `model_gc.json` and `model_1bin.json`;
- `model_gc_check.json`, which verifies that GC bins actually contain
  different transition structure;
- matched simulated FASTQ files for sample, one-bin, and GC-aware routes;
- `gc_improvement_summary.csv`, `.md`, and `.json`.

The gate records whether GC-aware improves the k-mer, GC, and composite
scores. A failed gate is retained as a negative result rather than
rewritten.

## Teacher-framework figures

The four-layer table and innovation figure suite can be rebuilt from the
result tree with:

```bash
bash reproduce.sh --framework-figures
```

The figure generator refuses to overwrite an existing output directory.
Plotted heatmap cells pass a strict finite-value gate. A tool is omitted from
a panel when the corresponding measurement is structurally unavailable, with
the exclusion stated in the caption, rather than being drawn as `NA`.

## Public data

The analyzed panel is defined in `config/species_config_v2.json`. ONT read
accessions and reference assembly accessions are recorded in
`config/data_sources.json`, `docs/methods_parameters.md`, and Table 1.
Large FASTQ/BAM/reference files are intentionally excluded from Git.

The six analyzed panels are:

| Species panel | ONT read accession |
|---|---|
| E. coli | SRR39619343 |
| S. cerevisiae | SRR39791672 |
| A. thaliana | ERR5716408 |
| D. melanogaster | SRR22071686 |
| M. musculus chr19 | SRR14685232 |
| H. sapiens chr21 | ERR13491992 |

## Full pipeline

The coverage-matched chr21 30x decision has its own provenance contract:

```bash
PROJECT_DIR=/path/to/project bash reproduce.sh --audit-panel
```

This independently checks that the decision, standalone figure, and exported
framework CSV rows all resolve to the `*_Hsapiens_chr21_30x.json` files. The
export step aborts rather than falling back to the broader `_alignment` or
`_fair_*` panel. See `docs/panel_provenance.md`.

Dependency-free checks:

```bash
bash reproduce.sh --tests
bash reproduce.sh --static-check
```

Optional release audits:

```bash
PROJECT_DIR=/path/to/project JOBS=4 bash reproduce.sh --audit-data
PROJECT_DIR=/path/to/project bash reproduce.sh --audit-logs
bash reproduce.sh --audit-paths
```

`--audit-data` checks BAM/gzip/FASTQ integrity and configured source hashes.
Known transient files and converted source FASTQs are reported separately
from critical failures. `--audit-logs` scans recent logs for failure signals,
and `--audit-paths` reports non-portable absolute paths.

The main implementation entry points are:

- `src/train_errhmm.py`
- `src/generate_simulated.py`
- `src/evaluate_level1.py`
- `src/downstream_tasks.py`
- `src/ablation_study.py`
- `scripts/cross_tool_compare.sh`

The full execution requires downloading the public data and constructing the
expected `data/` layout documented in `config/species_config_v2.json`.

For real HiFi validation, `config/hifi_sources.json` records ENA/SRA
accessions, official file sizes, and MD5 checksums. The download and
validation entry points are `scripts/titan_fetch_hifi.sh` and
`scripts/titan_validate_hifi.py`.

## Citation and DOI

The GitHub repository URL and Zenodo DOI must be filled after publication.
They are deliberately not fabricated here.
