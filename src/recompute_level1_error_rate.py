#!/usr/bin/env python3
"""Recompute Level-1 error-rate metrics for all existing Level-1 tables.

The pipeline originally disabled error-rate alignment because whole-FASTQ
minimap2 runs exhausted the 7.4GB WSL memory. This script samples the first
5000 reads per FASTQ before aligning, which keeps memory bounded while making
the error-rate subscore available for every species.
"""

import argparse
import json
import sys
import tempfile
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from evaluate_level1 import Level1Evaluator  # noqa: E402


ROUTES = ('route_A_sample', 'route_B_qshmm', 'route_C_errhmm')
CACHE_VERSION = 'mapq20_primary_reservoir_v3'


def parse_table_name(name):
    stem = Path(name).stem
    if stem.startswith('level1_'):
        stem = stem[len('level1_'):]
    for route in ROUTES:
        marker = f'_{route}_r'
        if marker in stem:
            species = stem.split(marker)[0]
            rep = stem.split(marker)[1]
            return species, route, int(rep)
    raise ValueError(f'Cannot parse Level-1 table name: {name}')


def real_fastq_for(project_dir, species):
    base = project_dir / 'data' / 'real_reads_verified'
    path = base / f'{species}_ont.fastq.gz'
    if not path.exists() and species == 'Hsapiens_chr21':
        path = base / 'Hsapiens_ont.fastq.gz'
    elif not path.exists() and species == 'Mmusculus_chr19':
        path = base / 'Mmusculus_ont.fastq.gz'
    return path


def load_json(path):
    with open(path, encoding='utf-8') as f:
        return json.load(f)


def save_json(path, data):
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2)


def compute_error_rate(evaluator, fastq, ref, threads, max_reads):
    with tempfile.TemporaryDirectory(prefix='level1_err_') as tmp_dir:
        return evaluator.align_and_compute_error_rate(
            str(fastq), str(ref), tmp_dir, threads=threads, max_reads=max_reads
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--project-dir', default=str(Path(__file__).resolve().parents[1]))
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--max-reads', type=int, default=5000)
    parser.add_argument('--min-mapq', type=int, default=20)
    parser.add_argument(
        '--species',
        default='',
        help='Comma-separated species subset, e.g. Ecoli,Scerevisiae',
    )
    args = parser.parse_args()

    project_dir = Path(args.project_dir)
    tables_dir = project_dir / 'results' / 'tables'
    cache_dir = project_dir / 'data' / 'tmp' / 'level1_error_rate_cache'
    cache_dir.mkdir(parents=True, exist_ok=True)

    evaluator = Level1Evaluator(min_mapq=args.min_mapq)
    table_files = sorted(tables_dir.glob('level1_*_r?.json'))
    species_filter = {
        item.strip() for item in args.species.split(',') if item.strip()
    }
    if species_filter:
        filtered = []
        for table_file in table_files:
            try:
                species, _, _ = parse_table_name(table_file.stem)
            except ValueError:
                continue
            if species in species_filter:
                filtered.append(table_file)
        table_files = filtered
    updated = []
    skipped = []

    real_cache = {}
    for table_file in table_files:
        try:
            species, route, rep = parse_table_name(table_file.stem)
        except ValueError as exc:
            print(f'[WARN] {exc}')
            skipped.append(str(table_file))
            continue

        ref = project_dir / 'data' / 'reference' / f'{species}_ref.fa'
        real_fastq = real_fastq_for(project_dir, species)
        sim_fastq = (
            project_dir / 'data' / 'simulated' / route / species /
            f'{species}_{route}_r{rep}.fastq.gz'
        )
        if not ref.exists() or not real_fastq.exists() or not sim_fastq.exists():
            print(f'[SKIP] missing inputs: {table_file.name}')
            skipped.append(str(table_file))
            continue

        results = load_json(table_file)
        old_composite = results.get('composite_score', {})
        if 'composite_score_4d' not in results:
            results['composite_score_4d'] = old_composite.get('composite_score')
        if 'sub_scores_4d' not in results:
            results['sub_scores_4d'] = old_composite.get('sub_scores')

        if species not in real_cache:
            source_stat = real_fastq.stat()
            cache_file = cache_dir / (
                f'{species}_real_{source_stat.st_size}_'
                f'{source_stat.st_mtime_ns}_{CACHE_VERSION}.json'
            )
            if cache_file.exists():
                real_cache[species] = load_json(cache_file)
                print(f'[CACHE] real error rate: {species}')
            else:
                real_cache[species] = compute_error_rate(
                    evaluator, real_fastq, ref, args.threads, args.max_reads
                )
                if real_cache[species] is not None:
                    real_cache[species]['source_size'] = source_stat.st_size
                    real_cache[species]['source_mtime_ns'] = source_stat.st_mtime_ns
                    save_json(cache_file, real_cache[species])
                    print(f'[DONE] real error rate: {species}')
                else:
                    print(f'[FAIL] real error rate: {species}')

        sim_error = compute_error_rate(
            evaluator, sim_fastq, ref, args.threads, args.max_reads
        )
        real_error = real_cache.get(species)
        if real_error is None or sim_error is None:
            print(f'[SKIP] error-rate failed: {table_file.name}')
            skipped.append(str(table_file))
            continue

        comparison = results.setdefault('comparison', {})
        denom = max(real_error['total'], sim_error['total'])
        comparison['error_rate_similarity'] = (
            max(0.0, 1.0 - abs(real_error['total'] - sim_error['total']) / denom)
            if denom > 0 else 0.0
        )
        comparison['real_error_rate'] = real_error['total']
        comparison['sim_error_rate'] = sim_error['total']
        comparison['real_error_rate_components'] = real_error
        comparison['sim_error_rate_components'] = sim_error
        comparison['error_rate_metric'] = CACHE_VERSION
        results['composite_score'] = evaluator.compute_composite_score(comparison)
        save_json(table_file, results)
        updated.append(table_file.name)
        print(f'[DONE] {table_file.name}')

    manifest = {
        'generated_at': datetime.now().isoformat(timespec='seconds'),
        'error_rate_mode': CACHE_VERSION,
        'error_rate_filter': {
            'primary_alignments_only': True,
            'min_mapq': args.min_mapq,
            'min_aligned_query_bases': 100,
        },
        'max_reads_per_file': args.max_reads,
        'threads': args.threads,
        'error_rate_weight': 0.25,
        'tables_updated': updated,
        'tables_skipped': skipped,
    }
    save_json(tables_dir / 'level1_error_rate_manifest.json', manifest)
    print(f'[SUMMARY] updated={len(updated)} skipped={len(skipped)}')


if __name__ == '__main__':
    main()
