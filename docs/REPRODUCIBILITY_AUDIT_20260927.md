# Reproducibility audit, 2026-09-27

## Scope

The audit covered:

* Python and shell syntax across `src/`, `scripts/`, and `tests/`;
* dependency-free unit tests;
* frozen artifact hashes in the latest block of each path in
  `docs/data_freeze_manifest.md`;
* the coverage-matched chr21 30x decision panel and its exported CSV
  derivatives;
* repeatability of the six main figures in PNG and PDF;
* BAM, gzip, FASTQ and configured-source integrity on Titan;
* recent pipeline logs;
* hardcoded-path and repository-cleanliness findings.

## Passed checks

| Check | Result |
|---|---|
| Python/shell syntax | `STATIC_AUDIT_OK python_files=92 shell_files=158` |
| Unit tests | `TESTS_OK passed=11` |
| Frozen artifacts | `FREEZE_MANIFEST checked=61 failed=0 missing=0` |
| 30x provenance | `PANEL_PROVENANCE_OK`, GCerrHMM R3 rank 2.0 and R5 rank 1.0 |
| Derived tables | Re-running export + decision + layer summary produced byte-identical outputs |
| Main figures | PNG and PDF hashes identical across two independent builds |
| GC improvement reproduction | Exit 0; GC-aware minus 1-bin composite +0.23, k-mer +1.08, GC +0.05 |
| Data integrity | 187 BAM + 529 gzip + 290 FASTQ files checked; 0 critical failures |

The deterministic main-figure snapshot is
`docs/gcerrhmm_main_figures_final_20260927_2145/`. Figure 1c now contains the
six-species reference/ONT table from `table1_species_panel.csv`; the
previous "Species table unavailable" branch is removed and missing input now
raises an error. The manuscript skeleton
`docs/manuscript/GCerrHMM_论文骨架_W1.3.1.docx` has SHA-256
`505c3639daf2c27c1ca88fa4fb5fa0b101692289a5a2dbafe8bc7234b1df988a`
and its key values match the 30x decision: R3 rank 2, R5 rank 1, top-2 in
2/2 anchored layers, chr21 assembly N50 25.04 Mbp against a 33.50 Mbp real
ONT anchor.

The end-to-end `--gc-demo` routes scored:

| route | composite | read length | QV | GC | k-mer |
|---|---:|---:|---:|---:|---:|
| sample_A | 85.31 | 99.24 | 70.09 | 97.90 | 70.72 |
| errhmm_1bin | 86.69 | 99.19 | 70.10 | 98.68 | 76.32 |
| errhmm_gc | 86.92 | 99.19 | 70.10 | 98.73 | 77.40 |

The reported gate was positive on composite, GC and k-mer.

## Data-integrity warnings

The audit deliberately separates critical failures from known non-input
intermediates. The final run reported:

```
checked=1020 critical_failures=0 known_warnings=9
AUDIT_EXIT=0
```

Warnings:

* one empty BAM under `data/tmp/`;
* one empty smoke-test VCF under `results/framework/phasing_ecoli_smoke/`;
* one empty placeholder `HG002_chr21.fastq`;
* two superseded FASTQ files under the old `cross_tools` panel;
* source-object MD5s for three ONT files differ from the recompressed local
  conversion; local gzip integrity passes;
* the optional human HiFi source `Hsapiens_ERR13110527_hifi.fastq.gz` was
  incomplete/corrupt; it has been moved with its `.aria2` control file to
  `data/real_reads_hifi/_incomplete_Hsapiens_ERR13110527_20260927/`.

The optional human HiFi source is not used by the default manuscript HiFi
chain, which uses `Ecoli_ERR12723508_hifi.fastq.gz`. A future human-HiFi
branch must re-download and validate it before use.

## Logs and repository hygiene

The three-day log scan covered 221 files and found 81 historical
failure-signal groups. The important CRC32 failures for A. thaliana and
M. musculus were followed by re-alignment; the final 48/48 GC-fidelity table
is complete, and all 187 BAMs passed `samtools quickcheck` in the final data
audit. The log scan is therefore diagnostic, not a substitute for the
artifact-level gates.

The hardcoded-path audit found 57 occurrences in 249 code/script files, mostly
in historical WSL-local drivers that are not part of the public
`reproduce.sh` path. The public path is parameterized through `PROJECT_DIR`.
These local-only scripts remain a release-cleanup item.

## Commands

The reproducible checks are now exposed as:

```bash
bash reproduce.sh --tests
bash reproduce.sh --static-check
PROJECT_DIR=/path/to/project bash reproduce.sh --full-check
PROJECT_DIR=/path/to/project bash reproduce.sh --audit-panel
PROJECT_DIR=/path/to/project bash reproduce.sh --audit-data
PROJECT_DIR=/path/to/project bash reproduce.sh --audit-logs
bash reproduce.sh --audit-paths
```

## Residual risk

The core paper path passes all executable gates. Remaining risks are:

1. the optional human HiFi source is archived as incomplete and would need a
   fresh download before use;
2. historical local-WSL driver scripts still contain hardcoded paths;
3. the repository contains old proposal/render duplicates that should be
   excluded or moved during public-release curation.
