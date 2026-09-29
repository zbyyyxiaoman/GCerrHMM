# GCerrHMM

**A reproducible GC-aware error HMM for Oxford Nanopore read simulation.**

GCerrHMM learns error-state transitions directly from a real ONT alignment and
conditions the transition matrices on local GC content. The generated reads
therefore preserve GC-dependent error structure, not only the global average
error rate.

This repository is the reviewer-facing implementation package. It contains the
model, the simulation and evaluation harness, the pre-registered decision
rules, the frozen result checks, and the manuscript/figure builders. It does
not redistribute sequencing reads.

## Why this repository exists

The central question is not whether a simulator can match one aggregate score.
It is whether a model changes the *right* quantity when the reference has
non-uniform GC structure.

The evaluation therefore separates two things:

* **GC-conditioned fidelity**: a GC-stratified error curve, scored by Pearson
  correlation and mean absolute deviation.
* **Overall descriptive fidelity**: a composite score over read length, QV,
  GC distribution and k-mer composition. This score is useful as a summary,
  but it does not by itself test the GC-conditioning claim.

This distinction is deliberate. A high composite score cannot substitute for
a GC-specific measurement.

## Scientific status

The current evidence is directional for a species-dependent effect:

* In *A. thaliana*, five or more GC bins improve the shape fidelity of the
  local GC-error curve over the one-bin control: `r = 0.779` at one bin versus
  `r >= 0.960` at five or more bins (`n = 3` seeds).
* The six-species panel favours the GC-aware route in four of six genomes in
  both replicates. The exceptions are *S. cerevisiae* (narrow GC structure)
  and *E. coli* (limited discriminating power at this coverage).
* The strict pre-registered two-species gate was not met. The result is
  therefore reported as directional evidence plus effect size, not as a
  universal claim that GC conditioning always improves every simulator
  metric.
* In the coverage-matched human chr21 30x panel, GCerrHMM ranks second in
  alignment identity (R3) and first in the assembly Borda comparison (R5),
  so it is top-two in both layers that have a structurally valid real anchor.
  R4 is excluded because its truth set is synthetic.

The exact decision rules, negative results and provenance are recorded under
[`docs/`](docs/). The manuscript package is under
[`docs/manuscript/`](docs/manuscript/).

For a step-by-step reviewer checklist, see
[`docs/REVIEWER_GUIDE.md`](docs/REVIEWER_GUIDE.md).

## Quick start

### 1. Create the core environment

```bash
conda env create -f environment.yml
conda activate errhmm
```

The external simulator comparison environment is separate:

```bash
conda env create -f environment-tools.yml
conda activate toolcompare
```

### 2. Run the data-free checks first

These commands work from a clean clone and do not require sequencing reads:

```bash
bash reproduce.sh --tests
bash reproduce.sh --static-check
bash reproduce.sh --smoke
```

`--smoke` verifies the SHA-256 checksums of the frozen result tables and
figure archives. `--tests` runs the dependency-light unit tests.
`--static-check` compiles the Python entry points and syntax-checks the shell
launcher.

### 3. Run the fastest end-to-end algorithmic check

The quick demo requires one real reference, one aligned BAM and the matching
FASTQ at the documented paths below. It trains the GC-aware model and the
one-bin control, simulates three matched routes, evaluates all three routes,
and writes the C-minus-B report:

```bash
bash reproduce.sh --quick-demo
```

Quick mode uses deterministic 1x coverage, reduced evaluation sampling and
parallel training/simulation/evaluation jobs. It is designed to finish in
minutes, not hours, while keeping the same code path as the full demo. It is a
smoke test, not the publication-scale result.

For the publication-scale 10x run:

```bash
SPECIES=Ecoli THREADS=16 JOBS=3 bash reproduce.sh --gc-demo
```

Both commands respect a 16-thread ceiling. The launcher applies
`nice -n 19` on systems where `nice` is available. Existing run directories are
never overwritten; each run gets a timestamped output directory.

### 4. Verify the frozen project tree when the full data are available

```bash
PROJECT_DIR=/path/to/project bash reproduce.sh --full-check
```

`--full-check` verifies every hash in
[`docs/data_freeze_manifest.md`](docs/data_freeze_manifest.md). It is expected
to fail in a bare clone because the large input and result trees are
intentionally not redistributed.

## Data layout

The code accepts an external project root through `PROJECT_DIR`. The demo
paths below are relative to that root:

```text
data/
  references/
    <species>_ref.fa
  truth/
    <species>_variant_ref.fa          # optional; falls back to the reference
  real_reads_verified/
    <species>_ont_aligned.bam
    <species>_ont.fastq.gz
  tmp/
    read_profiles/
      <species>.json                  # optional external profile
```

Read-length and QV profiles are deterministic metadata. If an external
profile is absent, the tracked fallback under
`config/read_profiles/<species>.json` is used. The fallback contains no reads
and no machine-specific paths; it makes the quick demo reproducible across
machines.

Public read accessions, expected sizes and checksums are recorded in
[`config/data_sources.json`](config/data_sources.json). The download helpers
refuse to continue on a checksum mismatch.

## Commands

| Command | Purpose |
|---|---|
| `bash reproduce.sh --smoke` | Verify frozen artifact checksums without rerunning science. |
| `bash reproduce.sh --tests` | Run the dependency-light unit tests. |
| `bash reproduce.sh --static-check` | Compile Python files and syntax-check the shell launcher. |
| `bash reproduce.sh --quick-demo` | Fast deterministic one-species train/simulate/evaluate chain. |
| `bash reproduce.sh --gc-demo` | Publication-scale 10x train/simulate/evaluate chain. |
| `bash reproduce.sh --r1` | Rebuild the six-species Level-1 summary from frozen JSON. |
| `bash reproduce.sh --audit-panel` | Verify that the chr21 30x decision and figures use the pre-registered panel. |
| `bash reproduce.sh --framework-figures` | Rebuild summary tables and framework figures. |
| `bash reproduce.sh --full-check` | Verify the full data freeze manifest when the project data are present. |
| `bash reproduce.sh --audit-paths` | Scan the release for non-portable absolute paths; set `PATH_AUDIT_STRICT=1` to fail on findings. |

## Outputs

The quick/full demo writes:

```text
results/gc_improvement/<species>_<timestamp>/
  model_gc.json
  model_1bin.json
  model_gc_check.json
  simulated/
  evaluation_<route>.json
  gc_improvement_summary.json
  gc_improvement_summary.csv
  gc_improvement_summary.md
  cache/
  run.log
  logs/
```

The summary records the C-minus-B deltas for composite, read-length, QV, GC and
k-mer metrics. A failed gate is preserved as a negative result; it is never
rewritten into a positive claim.

## Repository layout

```text
src/
  train_errhmm.py          HMM training from a real BAM
  generate_simulated.py   deterministic read simulation and coverage gate
  evaluate_level1.py      length/QV/GC/k-mer evaluation
  evaluate_sim_fastq.py   one-route Level-1 evaluation entry point
  downstream_tasks.py     mapping/variant/assembly/phasing helpers
scripts/
  reproduce_gc_improvement.py   end-to-end quick/full demo
  delta_to_real_decision.py     pre-registered chr21 decision rule
  build_*_figures.py            manuscript figure builders
  package_github_release.py     clean portable release builder
config/
  data_sources.json             public accessions and checksums
  read_profiles/                deterministic length/QV fallbacks
docs/
  delta_to_real_framing.md      pre-registered decision rule
  data_freeze_manifest.md       append-only artifact hashes
  panel_provenance.md           30x panel provenance contract
  manuscript/                   W3.3 submission and review copies
tests/
  run_tests.py                   dependency-free test entry point
```

## Reproducibility contract

1. **No hidden input fallback for measured runs.** Missing or short outputs
   fail loudly instead of silently changing the simulation regime.
2. **One decision source per claim.** The coverage-matched `*_30x.json` panel
   is authoritative for the delta-to-real decision.
3. **Frozen artifacts are append-only.** New analyses use a new timestamped
   output; existing hashes are not overwritten.
4. **Portability is checked.** The release builder rejects private usernames,
   private hostnames and machine-local paths.
5. **Resource use is bounded.** Long-running jobs use `nice -n 19` and a
   maximum of 16 worker threads across the project.

## Manuscript and citation

The current submission package is under
[`docs/manuscript/`](docs/manuscript/). Citation metadata is in
[`CITATION.cff`](CITATION.cff).

Authors:

* Boyang Zhang, Faculty of Life Science and Medicine, Harbin Institute of
  Technology.
* Tao Jiang, Faculty of Computing, Harbin Institute of Technology
  (corresponding author).

The repository URL is `https://github.com/zbyyyxiaoman/GCerrHMM`; the Zenodo
DOI is filled after archival release.

## Licence

MIT. See [`LICENSE`](LICENSE).
