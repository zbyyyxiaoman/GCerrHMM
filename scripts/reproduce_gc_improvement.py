#!/usr/bin/env python3
"""Reproduce the core GC-aware improvement loop on one species.

The script trains two models from the same real BAM, simulates matched-coverage
reads through three routes, evaluates each route against the same real FASTQ,
and writes a machine-readable improvement report. It never overwrites frozen
results.

Use ``--quick`` for a fast end-to-end smoke run with deterministic reduced
sampling. The full defaults are intended for the publication-scale check.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path


THREAD_CAP = 16

FULL_SETTINGS = {
    "coverage": 10.0,
    "max_reads": 10000,
    "evaluation_reads": None,
    "kmer_reads": None,
    "kmer_sample_size": None,
}

QUICK_SETTINGS = {
    "coverage": 1.0,
    "max_reads": 1000,
    "evaluation_reads": 500,
    "kmer_reads": 50,
    "kmer_sample_size": 2000,
}


def run(command, log_path, env=None):
    printable = " ".join(str(item) for item in command)
    print("[GC-IMPROVE] " + printable, flush=True)
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("a", encoding="utf-8") as log:
        log.write("$ " + printable + "\n")
        log.flush()
        result = subprocess.run(
            [str(item) for item in command],
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
            env=env,
        )
    if result.returncode != 0:
        raise RuntimeError(
            f"command failed ({result.returncode}): {command}; log={log_path}"
        )


def nice_prefix(enabled):
    return ["nice", "-n", "19"] if enabled and shutil.which("nice") else []


def workload_environment():
    """Prevent nested BLAS pools from exceeding the explicit thread budget."""
    env = os.environ.copy()
    for name in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        env[name] = "1"
    return env


def resolve_run_settings(
    quick: bool,
    coverage=None,
    max_reads=None,
    evaluation_reads=None,
    kmer_reads=None,
    kmer_sample_size=None,
):
    """Resolve full/quick defaults while letting explicit CLI values win."""
    defaults = QUICK_SETTINGS if quick else FULL_SETTINGS
    explicit = {
        "coverage": coverage,
        "max_reads": max_reads,
        "evaluation_reads": evaluation_reads,
        "kmer_reads": kmer_reads,
        "kmer_sample_size": kmer_sample_size,
    }
    return {
        key: defaults[key] if value is None else value
        for key, value in explicit.items()
    }


def resolve_threads(requested=None):
    """Return a positive thread budget capped at the shared-machine limit."""
    if requested is None:
        requested = min(THREAD_CAP, os.cpu_count() or 1)
    return max(1, min(THREAD_CAP, int(requested)))


def allocate_threads(total: int, workers: int) -> list[int]:
    """Split a thread budget across concurrent workers without oversubscribing."""
    total = max(1, int(total))
    workers = max(1, min(int(workers), total))
    base, extra = divmod(total, workers)
    return [
        base + (1 if index < extra else 0)
        for index in range(workers)
    ]


def run_commands(items, jobs, log_dir, env):
    """Run independent commands concurrently with one log file per command."""
    log_dir = Path(log_dir)

    def execute(item):
        label, command = item
        path = log_dir / f"{label}.log"
        run(command, path, env=env)
        return label, path

    if jobs <= 1 or len(items) <= 1:
        return [execute(item) for item in items]

    outcomes = []
    failures = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as executor:
        futures = [executor.submit(execute, item) for item in items]
        for future in concurrent.futures.as_completed(futures):
            try:
                outcomes.append(future.result())
            except Exception as exc:  # keep waiting for the other jobs
                failures.append(exc)
    if failures:
        raise failures[0]
    order = {label: index for index, (label, _) in enumerate(items)}
    return sorted(outcomes, key=lambda item: order[item[0]])


def append_logs(run_log: Path, task_logs) -> None:
    with run_log.open("a", encoding="utf-8") as aggregate:
        for label, path in task_logs:
            aggregate.write(f"\n## {label}\n\n")
            if path.exists():
                aggregate.write(path.read_text(encoding="utf-8"))
            aggregate.write("\n")


def load_metric(path):
    payload = json.loads(path.read_text(encoding="utf-8"))
    composite = payload.get("composite_score", {})
    return {
        "composite": float(composite.get("composite_score", 0.0)),
        "read_length": float(composite.get("sub_scores", {}).get("read_length", 0.0)),
        "qv": float(composite.get("sub_scores", {}).get("qv", 0.0)),
        "gc": float(composite.get("sub_scores", {}).get("gc", 0.0)),
        "kmer": float(composite.get("sub_scores", {}).get("kmer", 0.0)),
    }


def write_report(output_dir, rows, metadata):
    csv_path = output_dir / "gc_improvement_summary.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["route", "composite", "read_length", "qv", "gc", "kmer"],
        )
        writer.writeheader()
        writer.writerows(rows)

    by_route = {row["route"]: row for row in rows}
    deltas = {
        metric: by_route["errhmm_gc"][metric] - by_route["errhmm_1bin"][metric]
        for metric in ("composite", "read_length", "qv", "gc", "kmer")
    }
    gate = {
        "kmer_positive": deltas["kmer"] > 0,
        "composite_positive": deltas["composite"] > 0,
        "gc_positive": deltas["gc"] > 0,
    }
    payload = {
        "metadata": metadata,
        "routes": rows,
        "gc_minus_1bin": deltas,
        "gate": gate,
    }
    (output_dir / "gc_improvement_summary.json").write_text(
        json.dumps(payload, indent=2),
        encoding="utf-8",
    )

    lines = [
        "# GC-aware improvement reproduction",
        "",
        f"- Species: `{metadata['species']}`",
        f"- Coverage: `{metadata['coverage']}x`",
        f"- Seed: `{metadata['seed']}`",
        f"- Simulation reference: `{metadata['simulation_ref']}`",
        f"- Run profile: `{'quick' if metadata['quick'] else 'full'}`",
        "",
        "| Route | Composite | Read length | QV | GC | k-mer |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['route']} | {row['composite']:.2f} | "
            f"{row['read_length']:.2f} | {row['qv']:.2f} | "
            f"{row['gc']:.2f} | {row['kmer']:.2f} |"
        )
    lines += [
        "",
        "## GC-aware minus one-bin control",
        "",
        "| Metric | Delta |",
        "|---|---:|",
    ]
    for metric, value in deltas.items():
        lines.append(f"| {metric} | {value:+.2f} |")
    lines += [
        "",
        "The gate is reported, not forced. A failed gate is preserved as a",
        "negative result and must not be rewritten into a positive claim.",
    ]
    (output_dir / "gc_improvement_summary.md").write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )
    return payload


def resolve_path(value, default):
    return Path(value).resolve() if value else default.resolve()


def display_path(path: Path, root: Path) -> str:
    """Store portable paths in generated metadata whenever possible."""
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return "<external-data>/" + path.name


def evaluation_command(
    repo: Path,
    real_fastq: Path,
    simulated: Path,
    reference: Path,
    output: Path,
    settings,
    cache_dir: Path,
):
    command = [
        sys.executable,
        repo / "src/evaluate_sim_fastq.py",
        "--real", real_fastq,
        "--sim", simulated,
        "--ref", reference,
        "--output", output,
        "--cache-dir", cache_dir,
    ]
    if settings["evaluation_reads"] is not None:
        command += ["--max-reads", str(settings["evaluation_reads"])]
    if settings["kmer_reads"] is not None:
        command += ["--kmer-reads", str(settings["kmer_reads"])]
    if settings["kmer_sample_size"] is not None:
        command += ["--kmer-sample-size", str(settings["kmer_sample_size"])]
    return command


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-dir", default=".")
    parser.add_argument(
        "--repo-dir",
        default=None,
        help="Directory containing src/ and scripts/; defaults to project-dir.",
    )
    parser.add_argument("--species", default="Ecoli")
    parser.add_argument("--reference", default=None)
    parser.add_argument("--simulation-reference", default=None)
    parser.add_argument("--bam", default=None)
    parser.add_argument("--real-fastq", default=None)
    parser.add_argument("--profile-json", default=None)
    parser.add_argument("--output-dir", default=None)
    parser.add_argument("--coverage", type=float, default=None)
    parser.add_argument("--platform", choices=("ont", "hifi"), default="ont")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help=f"Total worker budget; defaults to min({THREAD_CAP}, CPU count).",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=None,
        help="Concurrent training/simulation jobs; defaults to at most three.",
    )
    parser.add_argument(
        "--sequential",
        action="store_true",
        help="Disable stage-level concurrency for debugging.",
    )
    parser.add_argument("--gc-bins", type=int, default=10)
    parser.add_argument("--max-reads", type=int, default=None)
    parser.add_argument(
        "--evaluation-reads",
        type=int,
        default=None,
        help="Reservoir size for Level-1 read statistics.",
    )
    parser.add_argument(
        "--kmer-reads",
        type=int,
        default=None,
        help="Maximum reads used for the k-mer spectrum.",
    )
    parser.add_argument(
        "--kmer-sample-size",
        type=int,
        default=None,
        help="Maximum k-mer vocabulary retained for correlation.",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Use deterministic 1x / reduced-sampling settings for a fast smoke run.",
    )
    parser.add_argument("--nice", action="store_true", default=True)
    parser.add_argument("--no-nice", dest="nice", action="store_false")
    args = parser.parse_args()

    settings = resolve_run_settings(
        args.quick,
        coverage=args.coverage,
        max_reads=args.max_reads,
        evaluation_reads=args.evaluation_reads,
        kmer_reads=args.kmer_reads,
        kmer_sample_size=args.kmer_sample_size,
    )
    threads = resolve_threads(args.threads)
    jobs = args.jobs if args.jobs is not None else min(3, threads)
    jobs = max(1, min(int(jobs), 3))
    if args.sequential:
        jobs = 1

    root = Path(args.project_dir).resolve()
    repo = Path(args.repo_dir).resolve() if args.repo_dir else root
    data = root / "data"
    reference = resolve_path(
        args.reference, data / "references" / f"{args.species}_ref.fa"
    )
    simulation_ref = resolve_path(
        args.simulation_reference,
        data / "truth" / f"{args.species}_variant_ref.fa",
    )
    bam = resolve_path(
        args.bam, data / "real_reads_verified" / f"{args.species}_ont_aligned.bam"
    )
    real_fastq = resolve_path(
        args.real_fastq,
        data / "real_reads_verified" / f"{args.species}_ont.fastq.gz",
    )
    profile_candidates = [
        data / "tmp" / "read_profiles" / f"{args.species}.json",
        repo / "config" / "read_profiles" / f"{args.species}.json",
    ]
    if args.profile_json:
        profile_json = Path(args.profile_json).resolve()
    else:
        profile_json = next(
            (path for path in profile_candidates if path.exists()),
            profile_candidates[0],
        )
    if not simulation_ref.exists():
        simulation_ref = reference
    for label, path in (
        ("reference", reference),
        ("simulation reference", simulation_ref),
        ("real BAM", bam),
        ("real FASTQ", real_fastq),
    ):
        if not path.exists():
            raise SystemExit(
                f"missing {label}: {path}\n"
                "See README.md / REPRODUCIBILITY.md for the expected data layout."
            )

    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_root = (
        Path(args.output_dir).resolve()
        if args.output_dir
        else root / "results" / "gc_improvement"
    )
    output_dir = output_root / f"{args.species}_{run_id}_seed{args.seed}"
    output_dir.mkdir(parents=True, exist_ok=False)
    run_log = output_dir / "run.log"
    run_log.write_text(
        "GC improvement reproduction\n"
        f"profile={'quick' if args.quick else 'full'}\n"
        f"threads={threads}\n"
        f"jobs={jobs}\n\n",
        encoding="utf-8",
    )
    task_log_dir = output_dir / "logs"
    task_log_dir.mkdir()
    env = workload_environment()
    prefix = nice_prefix(args.nice)

    model_gc = output_dir / "model_gc.json"
    model_1bin = output_dir / "model_1bin.json"

    training = [
        (
            "train_gc",
            prefix + [
                sys.executable,
                repo / "src/train_errhmm.py",
                "--bam", bam,
                "--ref", reference,
                "--output", model_gc,
                "--gc-bins", str(args.gc_bins),
                "--max-reads", str(settings["max_reads"]),
            ],
        ),
        (
            "train_1bin",
            prefix + [
                sys.executable,
                repo / "src/train_errhmm.py",
                "--bam", bam,
                "--ref", reference,
                "--output", model_1bin,
                "--nongc",
                "--max-reads", str(settings["max_reads"]),
            ],
        ),
    ]
    training_logs = run_commands(
        training,
        jobs=min(jobs, len(training)),
        log_dir=task_log_dir,
        env=env,
    )
    append_logs(run_log, training_logs)

    check_log = task_log_dir / "check_gc_model.log"
    run(
        prefix + [
            sys.executable,
            repo / "scripts/check_gc_model.py",
            "--model", model_gc,
            "--output", output_dir / "model_gc_check.json",
            "--strict",
        ],
        check_log,
        env=env,
    )
    append_logs(run_log, [("check_gc_model", check_log)])

    simulate_root = output_dir / "simulated"
    simulate_root.mkdir()
    route_specs = [
        ("sample_A", None),
        ("errhmm_1bin", model_1bin),
        ("errhmm_gc", model_gc),
    ]
    route_jobs = min(jobs, len(route_specs))
    route_threads = allocate_threads(threads, route_jobs)
    route_commands = []
    for (route, model), thread_count in zip(route_specs, route_threads):
        route_prefix = simulate_root / f"{args.species}_{route}"
        command = prefix + [
            sys.executable,
            repo / "src/generate_simulated.py",
            "--ref", simulation_ref,
            "--output-prefix", route_prefix,
            "--coverage", str(settings["coverage"]),
            "--platform", args.platform,
            "--seed", str(args.seed),
            "--threads", str(thread_count),
        ]
        if model is not None:
            command += ["--errhmm-model", model]
        if profile_json.exists():
            command += ["--profile-json", profile_json]
        route_commands.append((f"simulate_{route}", command))

    simulation_logs = run_commands(
        route_commands,
        jobs=route_jobs,
        log_dir=task_log_dir,
        env=env,
    )
    append_logs(run_log, simulation_logs)

    evaluation_commands = []
    evaluation_cache = output_dir / "cache"
    for route, _ in route_specs:
        route_prefix = simulate_root / f"{args.species}_{route}"
        simulated = Path(f"{route_prefix}.fastq")
        evaluation = output_dir / f"evaluation_{route}.json"
        evaluation_commands.append(
            (
                f"evaluate_{route}",
                prefix + evaluation_command(
                    repo,
                    real_fastq,
                    simulated,
                    reference,
                    evaluation,
                    settings,
                    evaluation_cache,
                ),
            )
        )
    evaluation_logs = run_commands(
        evaluation_commands,
        jobs=route_jobs,
        log_dir=task_log_dir,
        env=env,
    )
    append_logs(run_log, evaluation_logs)

    summaries = []
    for route, _ in route_specs:
        summaries.append(
            {
                "route": route,
                **load_metric(output_dir / f"evaluation_{route}.json"),
            }
        )

    payload = write_report(
        output_dir,
        summaries,
        {
            "species": args.species,
            "coverage": settings["coverage"],
            "platform": args.platform,
            "seed": args.seed,
            "reference": display_path(reference, root),
            "simulation_ref": display_path(simulation_ref, root),
            "bam": display_path(bam, root),
            "real_fastq": display_path(real_fastq, root),
            "profile_json": (
                display_path(profile_json, root)
                if profile_json.exists()
                else None
            ),
            "gc_bins": args.gc_bins,
            "max_reads": settings["max_reads"],
            "evaluation_reads": settings["evaluation_reads"],
            "kmer_reads": settings["kmer_reads"],
            "kmer_sample_size": settings["kmer_sample_size"],
            "repo_dir": display_path(repo, root),
            "threads": threads,
            "jobs": jobs,
            "quick": bool(args.quick),
        },
    )
    print(f"GC_IMPROVEMENT_OUTPUT={output_dir}")
    print(json.dumps(payload["gc_minus_1bin"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
