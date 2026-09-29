#!/usr/bin/env python3
"""Run a real-read phasing evaluation with whatsHap.

The evaluator is intentionally explicit about ploidy. If the sample has no
heterozygous calls, it reports that condition instead of inventing a phase
block or switch-error value.
"""

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path


def run(command, log, check=True):
    print("[PHASING] " + " ".join(str(item) for item in command), flush=True)
    log.write("$ " + " ".join(str(item) for item in command) + "\n")
    log.flush()
    result = subprocess.run(
        [str(item) for item in command],
        stdout=log,
        stderr=subprocess.STDOUT,
        text=True,
    )
    if check and result.returncode != 0:
        raise RuntimeError(f"command failed ({result.returncode}): {command}")
    return result


def nice_prefix():
    return ["nice", "-n", "19"] if shutil.which("nice") else []


def mapped_percent(bam, log):
    result = subprocess.run(
        ["samtools", "flagstat", str(bam)],
        capture_output=True, text=True, check=True,
    )
    log.write(result.stdout)
    for line in result.stdout.splitlines():
        if " mapped (" in line and "%" in line:
            return float(line.split("(")[1].split("%")[0])
    return 0.0


def vcf_count(vcf, expression=None, log=None):
    command = ["bcftools", "view", "-H"]
    if expression:
        command += ["-i", expression]
    command.append(str(vcf))
    result = subprocess.run(command, capture_output=True, text=True, check=True)
    if log is not None:
        log.write(" ".join(command) + f"\ncount={len(result.stdout.splitlines())}\n")
    return len(result.stdout.splitlines())


def count_heterozygous(vcf, log):
    result = subprocess.run(
        ["bcftools", "query", "-f", "[%GT]\n", str(vcf)],
        capture_output=True, text=True, check=True,
    )
    values = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    hetero = [
        value for value in values
        if re.search(r"(?:0[/|]1|1[/|]0)", value)
    ]
    log.write(f"heterozygous_calls={len(hetero)} total_genotypes={len(values)}\n")
    return len(hetero)


def parse_stats_tsv(path):
    if not path.exists():
        return {}
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line]
    if len(lines) < 2:
        return {}
    header = lines[0].split("\t")
    values = lines[1].split("\t")
    parsed = {}
    for key, value in zip(header, values):
        try:
            parsed[key] = float(value)
        except ValueError:
            parsed[key] = value
    return parsed


def parse_compare_tsv(path):
    if not path.exists():
        return {}
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line]
    if len(lines) < 2:
        return {}
    header = lines[0].split("\t")
    values = lines[1].split("\t")
    parsed = {}
    for key, value in zip(header, values):
        try:
            parsed[key] = float(value)
        except ValueError:
            parsed[key] = value
    return parsed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reads", required=True)
    parser.add_argument("--ref", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--platform", choices=("ont", "hifi"), default="hifi")
    parser.add_argument("--truth-vcf", default=None)
    parser.add_argument("--region", default=None)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    reads = Path(args.reads).resolve()
    reference = Path(args.ref).resolve()
    truth = Path(args.truth_vcf).resolve() if args.truth_vcf else None
    for label, path in (("reads", reads), ("reference", reference)):
        if not path.exists():
            raise SystemExit(f"missing {label}: {path}")
    if truth is not None and not truth.exists():
        raise SystemExit(f"missing truth VCF: {truth}")

    output_dir = Path(args.output_dir).resolve()
    if output_dir.exists() and any(output_dir.iterdir()):
        raise SystemExit(f"refusing to overwrite non-empty output dir: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "phasing.log"
    prefix = nice_prefix()

    with log_path.open("w", encoding="utf-8") as log:
        bam = output_dir / "aligned.bam"
        preset = "map-hifi" if args.platform == "hifi" else "map-ont"
        run(
            prefix + [
                "bash", "-o", "pipefail", "-c",
                f"minimap2 -ax {preset} -t {args.threads} "
                f"{reference} {reads} | "
                f"samtools sort -@ {min(args.threads, 4)} -m 1G -o {bam}",
            ],
            log,
        )
        run(["/bin/true"], log)
        subprocess.run(["samtools", "index", str(bam)], check=True)
        mapping_rate = mapped_percent(bam, log)

        calls = output_dir / "calls.vcf.gz"
        run(
            prefix + [
                "bash", "-o", "pipefail", "-c",
                f"bcftools mpileup -f {reference} -q 20 -Q 10 -Ou {bam} | "
                f"bcftools call -mv -Oz -o {calls}",
            ],
            log,
        )
        subprocess.run(["bcftools", "index", str(calls)], check=True)
        filtered = output_dir / "filtered.vcf.gz"
        run(
            prefix + [
                "bcftools", "view", "-i", "QUAL>=20 && INFO/DP>=4",
                "-Oz", "-o", filtered, calls,
            ],
            log,
        )
        subprocess.run(["bcftools", "index", str(filtered)], check=True)

        raw_variants = vcf_count(calls, log=log)
        filtered_variants = vcf_count(filtered, log=log)
        heterozygous = count_heterozygous(filtered, log)
        payload = {
            "sample": args.sample,
            "platform": args.platform,
            "reads": str(reads),
            "reference": str(reference),
            "mapping_rate": round(mapping_rate, 4),
            "raw_variants": raw_variants,
            "filtered_variants": filtered_variants,
            "heterozygous_calls": heterozygous,
            "status": "ready" if heterozygous > 0 else "no_heterozygous_variants",
        }

        if heterozygous > 0:
            phased = output_dir / "phased.vcf.gz"
            run(
                prefix + [
                    "whatshap", "phase",
                    "--reference", reference,
                    "--output", phased,
                    "--ignore-read-groups",
                    filtered, bam,
                ],
                log,
            )
            subprocess.run(["bcftools", "index", str(phased)], check=True)
            stats_tsv = output_dir / "phase_stats.tsv"
            run(
                prefix + [
                    "whatshap", "stats", "--tsv", stats_tsv, phased,
                ],
                log,
            )
            payload["phase_stats"] = parse_stats_tsv(stats_tsv)

            if truth is not None:
                truth_region = output_dir / "truth_region.vcf.gz"
                view_command = ["bcftools", "view"]
                if args.region:
                    view_command += ["-r", args.region]
                view_command += ["-Oz", "-o", truth_region, truth]
                run(prefix + view_command, log)
                subprocess.run(["bcftools", "index", str(truth_region)], check=True)
                truth_gt = subprocess.run(
                    ["bcftools", "query", "-f", "[%GT]\n", str(truth_region)],
                    capture_output=True, text=True, check=False,
                ).stdout
                if "|" not in truth_gt:
                    payload["comparison"] = {}
                    payload["comparison_status"] = "truth_unphased"
                else:
                    compare_tsv = output_dir / "compare.tsv"
                    run(
                        prefix + [
                            "whatshap", "compare",
                            "--ignore-sample-name",
                            "--only-snvs",
                            "--tsv-pairwise", compare_tsv,
                            truth_region, phased,
                        ],
                        log,
                        check=False,
                    )
                    payload["comparison"] = parse_compare_tsv(compare_tsv)
                    payload["comparison_status"] = (
                        "completed" if compare_tsv.exists() else "not_available"
                    )

        (output_dir / "phasing_result.json").write_text(
            json.dumps(payload, indent=2), encoding="utf-8"
        )
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
