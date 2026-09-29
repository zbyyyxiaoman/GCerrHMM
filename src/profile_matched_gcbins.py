#!/usr/bin/env python3
"""Run one GC-bin ablation point with a real read-length/QV profile."""

import argparse
import json
from pathlib import Path

from evaluate_level1 import compare_real_vs_simulated
from generate_simulated import GCAwareSimulator


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-dir", required=True)
    parser.add_argument("--species", required=True)
    parser.add_argument("--gc-bins", type=int, required=True)
    parser.add_argument("--model", required=True)
    parser.add_argument("--profile-json", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--run-tag", required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--threads", type=int, default=1,
                        help="Worker count for the simulation step")
    args = parser.parse_args()

    project = Path(args.project_dir)
    data_dir = project / "data"
    ref_fasta = data_dir / "references" / f"{args.species}_ref.fa"
    real_fastq = data_dir / "real_reads_verified" / f"{args.species}_ont.fastq.gz"
    output = Path(args.output)
    sim_prefix = (
        data_dir / "simulated" / "ablation"
        / f"{args.species}_{args.gc_bins}bin_profile_{args.run_tag}_seed{args.seed}"
    )

    for path in (ref_fasta, real_fastq, Path(args.model), Path(args.profile_json)):
        if not path.exists():
            raise FileNotFoundError(str(path))
    if output.exists():
        raise FileExistsError(str(output))

    simulator = GCAwareSimulator(
        errhmm_model=args.model,
        platform="ont",
        profile_json=args.profile_json,
        seed=args.seed,
    )
    sim_fastq = simulator.simulate_reads(
        str(ref_fasta),
        coverage=10,
        output_prefix=str(sim_prefix),
        threads=args.threads,
    )
    evaluation = compare_real_vs_simulated(
        str(real_fastq),
        sim_fastq,
        str(ref_fasta),
    )
    composite = evaluation.get("composite_score", {})
    output.write_text(json.dumps({
        "species": args.species,
        "gc_bins": args.gc_bins,
        "profile_json": args.profile_json,
        "seed": args.seed,
        "model": args.model,
        "composite_score": composite.get("composite_score", 0),
        "sub_scores": composite.get("sub_scores", {}),
    }, indent=2))
    print(
        f"PROFILE_GCBINS_DONE species={args.species} "
        f"gc_bins={args.gc_bins} score={composite.get('composite_score', 0)}"
    )


if __name__ == "__main__":
    main()
