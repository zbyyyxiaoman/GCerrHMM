#!/usr/bin/env python3
"""Add the gc_bins=1 route-B control to an existing gc-bin ablation JSON."""

import argparse
import json
import os
from pathlib import Path

from evaluate_level1 import compare_real_vs_simulated
from generate_simulated import GCAwareSimulator


def evaluate_species(project_dir, species, source_json, output_json, run_tag):
    project = Path(project_dir)
    data_dir = project / "data"
    ref_fasta = data_dir / "references" / f"{species}_ref.fa"
    real_fastq = data_dir / "real_reads_verified" / f"{species}_ont.fastq.gz"
    model_path = data_dir / "trained_models" / f"{species}_errhmm_nongc.json"
    sim_prefix = data_dir / "simulated" / "ablation" / f"{species}_1bin_{run_tag}"

    for path in (ref_fasta, real_fastq, model_path, source_json):
        if not path.exists():
            raise FileNotFoundError(str(path))
    if output_json.exists():
        raise FileExistsError(str(output_json))

    simulator = GCAwareSimulator(errhmm_model=str(model_path), platform="ont")
    sim_fastq = simulator.simulate_reads(
        str(ref_fasta),
        coverage=10,
        output_prefix=str(sim_prefix),
        threads=1,
    )
    evaluation = compare_real_vs_simulated(
        str(real_fastq),
        sim_fastq,
        str(ref_fasta),
    )
    composite = evaluation.get("composite_score", {})
    results = json.loads(source_json.read_text())
    results["1bin"] = {
        "composite_score": composite.get("composite_score", 0),
        "sub_scores": composite.get("sub_scores", {}),
    }
    scores = {
        key: float(results.get(key, {}).get("composite_score", 0) or 0)
        for key in ("1bin", "5bin", "10bin", "20bin")
    }
    results["best_config"] = max(scores, key=scores.get)
    results["conclusion"] = (
        f"{results['best_config']} provides the best balance of "
        "accuracy and granularity"
    )
    output_json.write_text(json.dumps(results, indent=2))
    print(
        f"GCBIN1_DONE species={species} "
        f"composite_score={results['1bin']['composite_score']} "
        f"best_config={results['best_config']} output={output_json}"
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--source-json", required=True)
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--run-tag", required=True)
    args = parser.parse_args()
    os.environ.setdefault("TRAIN_MAX_READS", "20000")
    evaluate_species(
        args.project_dir,
        args.species,
        Path(args.source_json),
        Path(args.output_json),
        args.run_tag,
    )


if __name__ == "__main__":
    main()
