#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="${PROJECT_DIR:-$ROOT}"
MODE="${1:---smoke}"

run_gc_demo() {
    local profile="${1:-full}"
    local -a gc_args
    gc_args=(
        --project-dir "$PROJECT_DIR"
        --repo-dir "${CODE_DIR:-$ROOT}"
        --species "${SPECIES:-Ecoli}"
        --seed "${SEED:-42}"
        --threads "${THREADS:-16}"
        --jobs "${JOBS:-3}"
        --output-dir "${GC_OUTPUT_DIR:-$PROJECT_DIR/results/gc_improvement}"
    )
    if [[ "$profile" == "quick" ]]; then
        gc_args+=(--quick)
    fi
    if [[ -n "${COVERAGE:-}" ]]; then
        gc_args+=(--coverage "$COVERAGE")
    fi
    if [[ -n "${MAX_READS:-}" ]]; then
        gc_args+=(--max-reads "$MAX_READS")
    fi
    python3 "$ROOT/scripts/reproduce_gc_improvement.py" "${gc_args[@]}"
}

prepare_demo_data() {
    local -a prepare_args
    prepare_args=(
        --project-dir "$PROJECT_DIR" \
        --species "${SPECIES:-Ecoli}" \
        --threads "${THREADS:-16}" \
        --parts "${DOWNLOAD_PARTS:-8}" \
        --chunk-size "${DOWNLOAD_CHUNK_SIZE:-16777216}" \
        --sra-threads "${SRA_THREADS:-4}"
        --sort-memory "${SORT_MEMORY:-1G}"
    )
    if [[ -n "${SORT_THREADS:-}" ]]; then
        prepare_args+=(--sort-threads "$SORT_THREADS")
    fi
    if [[ -n "${MIN_FREE_GB:-}" ]]; then
        prepare_args+=(--min-free-gb "$MIN_FREE_GB")
    fi
    python3 "$ROOT/scripts/prepare_demo_data.py" "${prepare_args[@]}"
}

run_fixture_demo() {
    python3 "$ROOT/scripts/reproduce_gc_improvement.py" \
        --project-dir "$ROOT/tests/fixtures/tiny_demo" \
        --repo-dir "$ROOT" \
        --species Ecoli \
        --profile-json "$ROOT/tests/fixtures/tiny_demo/data/tmp/read_profiles/Ecoli.json" \
        --coverage "${FIXTURE_COVERAGE:-1}" \
        --seed "${SEED:-42}" \
        --threads "${THREADS:-16}" \
        --jobs "${JOBS:-3}" \
        --output-dir "${FIXTURE_OUTPUT_DIR:-$ROOT/results/fixture_demo}"
}

check_bundled_results() {
    python3 "$ROOT/scripts/install_bundled_results.py" \
        --project-dir "$PROJECT_DIR" \
        --bundle-dir "$ROOT/docs/reproducibility/results_bundle" \
        --check-only
}

ensure_bundled_results() {
    python3 "$ROOT/scripts/install_bundled_results.py" \
        --project-dir "$PROJECT_DIR" \
        --bundle-dir "$ROOT/docs/reproducibility/results_bundle" \
        --quiet
}

case "$MODE" in
    --smoke)
        check_bundled_results
        if command -v sha256sum >/dev/null 2>&1; then
            (cd "$ROOT" && sha256sum -c docs/reproducibility/checksums.sha256)
        else
            (cd "$ROOT" && shasum -a 256 -c docs/reproducibility/checksums.sha256)
        fi
        echo "REPRODUCE_SMOKE_OK"
        ;;
    --r1)
        ensure_bundled_results
        python3 "$ROOT/scripts/build_r1_level1_overall.py" --project-dir "$PROJECT_DIR"
        ;;
    --full-check)
        if [[ ! -d "$PROJECT_DIR/data" && ! -d "$PROJECT_DIR/results" ]]; then
            echo "FULL_CHECK_REQUIRES_PROJECT_DATA: set PROJECT_DIR to a checkout with data/ and results/" >&2
            exit 2
        fi
        python3 "$ROOT/scripts/verify_freeze_manifest.py" --project-dir "$PROJECT_DIR" --strict
        ;;
    --tests)
        python3 "$ROOT/tests/run_tests.py"
        ;;
    --static-check)
        python3 "$ROOT/scripts/static_audit.py" --repo-root "$ROOT"
        ;;
    --gc-demo)
        run_gc_demo full
        ;;
    --quick-demo)
        run_gc_demo quick
        ;;
    --prepare-demo)
        prepare_demo_data
        ;;
    --quick-reproduce)
        prepare_demo_data
        run_gc_demo quick
        ;;
    --fixture-demo)
        run_fixture_demo
        ;;
    --framework-figures)
        ensure_bundled_results
        python3 "$ROOT/scripts/export_framework_tables.py" \
            --project-dir "$PROJECT_DIR" \
            --output-dir "${OUTPUT_DIR:-$PROJECT_DIR/results/framework/stats}"
        python3 "$ROOT/scripts/build_gcerrhmm_main_figures.py" \
            --project-dir "$PROJECT_DIR" \
            --output-dir "${MAIN_FIGURE_DIR:-$PROJECT_DIR/docs/gcerrhmm_main_figures}"
        python3 "$ROOT/scripts/build_teacher_framework_figures.py" \
            --project-dir "$PROJECT_DIR" \
            --output-dir "${FIGURE_DIR:-$PROJECT_DIR/docs/framework_figures}"
        python3 "$ROOT/scripts/audit_panel_provenance.py" \
            --project-dir "$PROJECT_DIR"
        ;;
    --delta-to-real)
        ensure_bundled_results
        python3 "$ROOT/scripts/delta_to_real_decision.py" \
            --project-dir "$PROJECT_DIR" \
            --exclude-layer R4
        python3 "$ROOT/scripts/build_delta_to_real_figure.py" \
            --project-dir "$PROJECT_DIR"
        python3 "$ROOT/scripts/audit_panel_provenance.py" \
            --project-dir "$PROJECT_DIR"
        ;;
    --audit-panel)
        ensure_bundled_results
        python3 "$ROOT/scripts/audit_panel_provenance.py" \
            --project-dir "$PROJECT_DIR"
        ;;
    --audit-data)
        python3 "$ROOT/scripts/audit_data_integrity.py" \
            --project-dir "$PROJECT_DIR" \
            --jobs "${JOBS:-4}" \
            --output "${DATA_AUDIT_OUTPUT:-$PROJECT_DIR/logs/data_integrity_audit.md}"
        ;;
    --audit-logs)
        python3 "$ROOT/scripts/audit_logs.py" \
            --project-dir "$PROJECT_DIR" \
            --output "${LOG_AUDIT_OUTPUT:-$PROJECT_DIR/logs/log_audit.md}" \
            --days "${LOG_AUDIT_DAYS:-3}"
        ;;
    --audit-paths)
        path_args=(
            --project-dir "$ROOT"
            --output "${PATH_AUDIT_OUTPUT:-$PROJECT_DIR/logs/hardcoded_path_audit.md}"
        )
        if [[ "${PATH_AUDIT_STRICT:-0}" == "1" ]]; then
            path_args+=(--strict)
        fi
        python3 "$ROOT/scripts/audit_hardcoded_paths.py" \
            "${path_args[@]}"
        ;;
    --install-results)
        ensure_bundled_results
        ;;
    --additional-files)
        ensure_bundled_results
        python3 "$ROOT/scripts/build_additional_files.py" \
            --project-dir "$PROJECT_DIR" \
            --package-dir "$ROOT" \
            --output-dir "${ADDITIONAL_OUTPUT_DIR:-$PROJECT_DIR/docs/additional_files}"
        ;;
    --help|-h)
        cat <<'EOF'
Usage: bash reproduce.sh [mode]

  --smoke       Verify the frozen publication artifact checksums.
  --r1          Rebuild the R1 six-species x three-route Level-1 table.
  --full-check  Verify every checksum recorded in data_freeze_manifest.md.
  --tests       Run dependency-free unit tests.
  --static-check
                Compile every Python file and syntax-check every shell entry.
  --gc-demo     Retrain GC-aware and 1-bin models, resimulate, and compare.
  --quick-demo  Same chain with deterministic 1x/reduced-sampling presets.
  --prepare-demo
                Download one species reference/FASTQ and build its aligned BAM.
  --quick-reproduce
                Prepare missing demo data, then run --quick-demo.
  --fixture-demo
                Run the tiny offline fixture without downloading public data.
  --framework-figures
                Rebuild the teacher-framework summary tables and figures.
  --delta-to-real
                Apply the pre-registered delta-to-real rule to the chr21 30x
                panel and redraw the anchor figure.
  --audit-panel
                Verify that the decision, figure table, and exported framework
                CSVs all resolve to the coverage-matched 30x JSON panel.
  --audit-data
                Check BAM/gzip/FASTQ integrity and configured source MD5s.
  --audit-logs
                Summarise recent errors and integrity signals in logs.
  --audit-paths
                Report non-portable absolute paths in code and scripts.
  --install-results
                Materialise the frozen public result bundle under results/.
  --additional-files
                Build Additional files 1-7 from the frozen result bundle.
EOF
        ;;
    *)
        echo "Unknown mode: $MODE" >&2
        exit 2
        ;;
esac
