# Reviewer reproduction guide

This guide gives the shortest reliable path through the GCerrHMM submission.
It separates three checks that answer different questions.

## Tier 0: inspect the frozen result package

No sequencing reads are needed.

```bash
conda env create -f environment.yml
conda activate errhmm

python tests/run_tests.py
bash reproduce.sh --static-check
bash reproduce.sh --smoke
```

Expected final signals:

```text
TESTS_OK
STATIC_AUDIT_OK
REPRODUCE_SMOKE_OK
```

`--smoke` verifies the frozen Level-1 summary, provenance table and figure
archives against `docs/reproducibility/checksums.sha256`.

## Tier 1: run the algorithmic smoke chain

This is the fastest end-to-end test of the method itself. It retrains the
GC-aware model and its one-bin control from one real BAM, simulates three
matched routes, and evaluates all three with the Level-1 metrics.

The one-command path downloads the reference and public FASTQ, checks the
configured checksum, builds the aligned BAM with minimap2/samtools, and runs
the same deterministic quick chain:

```bash
PROJECT_DIR=/path/to/project \
CODE_DIR=/path/to/gcerrhmm/repository \
THREADS=16 \
JOBS=3 \
bash /path/to/gcerrhmm/repository/reproduce.sh --quick-reproduce
```

To prepare the inputs separately:

```bash
SPECIES=Ecoli THREADS=16 bash reproduce.sh --prepare-demo
bash reproduce.sh --quick-demo
```

Quick mode uses 1x coverage, 1,000 training reads, a 500-read evaluation
reservoir and deterministic profile fallbacks from
`config/read_profiles/`. It is designed for correctness and turnaround time,
not for reproducing the publication-scale effect size.

The command prints a path beginning with:

```text
GC_IMPROVEMENT_OUTPUT=...
```

Inspect:

```text
gc_improvement_summary.md
gc_improvement_summary.json
model_gc_check.json
run.log
logs/
```

The important reproducibility property is that the three routes share the
same seed and coverage. A negative delta is a valid result and must remain
negative.

## Tier 2: reproduce the publication-scale decision

This requires the downloaded project tree and the frozen result tables. Use
the same `PROJECT_DIR` from Tier 1.

```bash
PROJECT_DIR=/path/to/project bash reproduce.sh --r1
PROJECT_DIR=/path/to/project bash reproduce.sh --audit-panel
PROJECT_DIR=/path/to/project bash reproduce.sh --framework-figures
PROJECT_DIR=/path/to/project bash reproduce.sh --full-check
```

Expected provenance checks:

```text
PANEL_PROVENANCE_OK
FREEZE_MANIFEST checked=... failed=0 missing=0
```

`--audit-panel` is the important scientific gate. It verifies that the
delta-to-real decision, the standalone figure input, and the exported
framework tables all resolve to the pre-registered Hsapiens chr21 30x panel.

## Reading the claims

| Question | Evidence | What it supports |
|---|---|---|
| Does GC conditioning change local error structure? | GC-stratified error curve (`r`, MAD) | Directional, species-dependent evidence |
| Does the model improve the aggregate composite score? | Level-1 composite | Descriptive fidelity only |
| Does it match the real downstream anchor? | R3/R5 delta-to-real ranks | Top-two consistency in 2/2 valid anchored layers |
| Is the effect universal across species? | Pre-registered two-species gate | No; the gate was not met |

The repository deliberately does not convert a failed universal gate into a
success claim. The supported statement is narrower: GC conditioning produces a
reproducible directional effect in GC-heterogeneous genomes, with the clearest
effect in *A. thaliana*.

## Common failure modes

**Missing reference/BAM/FASTQ.** The demo exits with the expected path. Set
`PROJECT_DIR` to the directory containing `data/`, or pass explicit paths to
`reproduce_gc_improvement.py`.

**Bare clone fails `--full-check`.** That is expected: the full manifest points
to large files that are intentionally not redistributed. Run `--smoke`,
`--tests` and `--static-check` in a bare clone.

**A run directory already exists.** Run IDs are timestamped; the script
refuses to overwrite an existing run. Choose a new `GC_OUTPUT_DIR` or wait one
second before retrying.

**The quick demo has a negative GC delta.** Quick mode is deliberately
resource-limited. It checks that the chain runs and that the gate is reported
honestly. The publication-scale comparison is Tier 2.

## Release sanity check

Before uploading or handing the repository to a reviewer:

```bash
python scripts/package_github_release.py \
  --source . \
  --output dist/gcerrhmm_github_20260929 \
  --zip dist/GCerrHMM_github_20260929.zip

PATH_AUDIT_STRICT=1 bash dist/gcerrhmm_github_20260929/reproduce.sh --audit-paths
```

The packager starts from a clean staging directory and fails on private
usernames, private hostnames or machine-local paths.
