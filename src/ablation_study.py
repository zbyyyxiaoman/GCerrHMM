#!/usr/bin/env python3
"""
消融实验模块
功能: 系统性地验证errHMM各组件的贡献

消融实验设计:
1. GC-bin粒度: 5-bin vs 10-bin vs 20-bin → 验证GC分箱策略
2. 训练数据量: 1x, 5x, 10x, 30x → 验证数据需求
3. 跨平台迁移: HiFi训练→ONT测试 → 验证泛化能力
4. 状态空间: 简化版(仅M/S/I/D) vs 完整版(含D1-D4+) → 验证状态设计
"""

import sys
import os
import json
import numpy as np
import argparse
from pathlib import Path
from Bio import SeqIO
from train_errhmm import GCAwareErrHMM
from generate_simulated import GCAwareSimulator
from evaluate_level1 import Level1Evaluator, compare_real_vs_simulated


class AblationStudy:
    """消融实验执行器"""
    
    def __init__(self, project_dir, threads=16):
        self.project_dir = Path(project_dir)
        self.data_dir = self.project_dir / 'data'
        self.src_dir = self.project_dir / 'src'
        self.results_dir = self.project_dir / 'results'
        self.threads = threads
        self.train_max_reads = int(os.environ.get('TRAIN_MAX_READS', '10000'))
        self.ablation_reads_per_x = int(os.environ.get('ABLATION_READS_PER_X', '2000'))
        self.seed = int(os.environ.get('ABLATION_SEED', '42'))
        self.run_tag = os.environ.get('ABLATION_RUN_TAG', '').strip()
        self.run_suffix = f'_{self.run_tag}' if self.run_tag else ''

    def _is_current_model(self, model_path):
        try:
            with open(model_path, encoding='utf-8') as fh:
                return int(json.load(fh).get('model_version', 1)) >= 2
        except Exception:
            return False
    
    def ablation_gc_bins(self, test_species='Ecoli', output_json=None):
        """
        消融实验1: GC-bin粒度
        对比 5-bin, 10-bin, 20-bin 的效果
        """
        print("=" * 60)
        print("[消融实验1] GC-bin粒度对比")
        print("=" * 60)
        
        ref_fasta = self.data_dir / f'references/{test_species}_ref.fa'
        real_fastq = self.data_dir / f'real_reads_verified/{test_species}_ont.fastq.gz'
        bam_file = self.data_dir / f'real_reads_verified/{test_species}_ont_aligned.bam'
        
        if not ref_fasta.exists():
            print(f"  ⚠ 参考基因组不存在: {ref_fasta}")
            return None
        
        if not bam_file.exists():
            print(f"  ⚠ BAM文件不存在: {bam_file}")
            return None
        
        results = {}
        
        for gc_bins in [5, 10, 20]:
            print(f"\n--- {gc_bins}-bin ---")
            
            # 训练errHMM
            model_path = self.data_dir / f'trained_models/{test_species}_errhmm_{gc_bins}bin{self.run_suffix}.json'
            
            hmm = GCAwareErrHMM(gc_bins=gc_bins, window_size=100)
            hmm.train_from_bam(
                str(bam_file),
                str(ref_fasta),
                min_mapq=20,
                max_reads=self.train_max_reads,
            )
            hmm.save(str(model_path))
            
            # 生成模拟数据
            sim_prefix = self.data_dir / f'simulated/ablation/{test_species}_{gc_bins}bin{self.run_suffix}'
            sim_prefix.parent.mkdir(parents=True, exist_ok=True)
            
            simulator = GCAwareSimulator(errhmm_model=str(model_path), seed=self.seed)
            sim_fastq = simulator.simulate_reads(
                str(ref_fasta), coverage=10, 
                output_prefix=str(sim_prefix), threads=self.threads
            )
            
            # Level-1评估
            if real_fastq.exists():
                eval_result = compare_real_vs_simulated(
                    str(real_fastq), sim_fastq, str(ref_fasta)
                )
                results[f'{gc_bins}bin'] = {
                    'composite_score': eval_result.get('composite_score', {}).get('composite_score', 0),
                    'sub_scores': eval_result.get('composite_score', {}).get('sub_scores', {})
                }
            else:
                results[f'{gc_bins}bin'] = {'note': 'real data not available'}
        
        # 找出最优配置
        scores = {k: v.get('composite_score', 0) for k, v in results.items()}
        best = max(scores, key=scores.get) if scores else None
        results['best_config'] = best
        results['conclusion'] = f'{best} provides the best balance of accuracy and granularity'
        
        if output_json:
            with open(output_json, 'w') as f:
                json.dump(results, f, indent=2)
        
        print(f"\n✓ GC-bin消融完成")
        print(f"  最优配置: {best}")
        return results
    
    def ablation_coverage(self, test_species='Ecoli', output_json=None):
        """
        消融实验2: 训练数据覆盖度
        对比 1x, 5x, 10x, 30x 训练数据的效果
        """
        print("=" * 60)
        print("[消融实验2] 训练数据覆盖度对比")
        print("=" * 60)
        
        ref_fasta = self.data_dir / f'references/{test_species}_ref.fa'
        real_fastq = self.data_dir / f'real_reads_verified/{test_species}_ont.fastq.gz'
        full_bam = self.data_dir / f'real_reads_verified/{test_species}_ont_aligned.bam'
        
        if not full_bam.exists():
            print(f"  ⚠ BAM文件不存在")
            return None
        
        results = {}
        
        genome_length = sum(len(record.seq) for record in SeqIO.parse(str(ref_fasta), 'fasta'))
        source_config = self.project_dir / 'config' / 'data_sources.json'
        source_base_count = 0
        if source_config.exists():
            with open(source_config, encoding='utf-8') as handle:
                source_data = json.load(handle)
            source_base_count = source_data['species'].get(test_species, {}).get('base_count', 0)
        raw_coverage = source_base_count / genome_length if genome_length else 60.0

        for cov in [1, 5, 10, 30]:
            print(f"\n--- {cov}x coverage ---")
            
            # 子采样BAM到指定覆盖度
            subset_bam = self.data_dir / f'tmp/{test_species}_{cov}x_subset{self.run_suffix}.bam'
            subset_bam.parent.mkdir(parents=True, exist_ok=True)
            
            # 使用samtools view -s 进行子采样
            fraction = cov / 60  # 假设原始60x
            fraction = min(1.0, cov / max(raw_coverage, 1.0))
            import subprocess
            fraction_text = f"{fraction:.8f}".split(".", 1)[1]
            subprocess.run(
                f'samtools view -s {self.seed}.{fraction_text} -b {full_bam} '
                f'> {subset_bam} && samtools index {subset_bam}',
                shell=True
            )
            
            # 训练errHMM
            model_path = self.data_dir / f'trained_models/{test_species}_errhmm_{cov}x{self.run_suffix}.json'
            
            hmm = GCAwareErrHMM(gc_bins=10, window_size=100)
            hmm.train_from_bam(
                str(subset_bam),
                str(ref_fasta),
                min_mapq=20,
                max_reads=self.ablation_reads_per_x * cov,
            )
            hmm.save(str(model_path))
            
            # 生成模拟数据并评估
            sim_prefix = self.data_dir / f'simulated/ablation/{test_species}_{cov}x{self.run_suffix}'
            sim_prefix.parent.mkdir(parents=True, exist_ok=True)
            
            simulator = GCAwareSimulator(errhmm_model=str(model_path), seed=self.seed)
            sim_fastq = simulator.simulate_reads(
                str(ref_fasta), coverage=10,
                output_prefix=str(sim_prefix), threads=self.threads
            )
            
            if real_fastq.exists():
                eval_result = compare_real_vs_simulated(
                    str(real_fastq), sim_fastq, str(ref_fasta)
                )
                results[f'{cov}x'] = {
                    'composite_score': eval_result.get('composite_score', {}).get('composite_score', 0),
                    'convergence': eval_result.get('composite_score', {}).get('sub_scores', {})
                }
            else:
                results[f'{cov}x'] = {'note': 'real data not available'}
        
        # 分析收敛性
        scores = {k: v.get('composite_score', 0) for k, v in results.items()}
        if scores:
            sorted_covs = sorted(scores.keys(), key=lambda x: int(x[:-1]))
            score_values = [scores[c] for c in sorted_covs]
            
            # 计算边际收益
            marginal_gains = [score_values[i+1] - score_values[i] 
                            for i in range(len(score_values)-1)]
            
            results['marginal_gains'] = dict(zip(
                [f'{sorted_covs[i]}→{sorted_covs[i+1]}' for i in range(len(sorted_covs)-1)],
                marginal_gains
            ))
            
            # 推荐配置 (边际收益 < 1%的点)
            optimal = None
            for i, gain in enumerate(marginal_gains):
                if gain < 1.0:
                    optimal = sorted_covs[i]
                    break
            if not optimal:
                optimal = sorted_covs[-1]
            
            results['recommended_coverage'] = optimal
            results['conclusion'] = f'{optimal} coverage provides the best cost-performance trade-off'
        
        if output_json:
            with open(output_json, 'w') as f:
                json.dump(results, f, indent=2)
        
        print(f"\n✓ 覆盖度消融完成")
        if 'recommended_coverage' in results:
            print(f"  推荐覆盖度: {results['recommended_coverage']}")
        return results
    
    def ablation_cross_platform(self, test_species='Ecoli', output_json=None):
        """
        消融实验3: 跨平台迁移
        按文档设计: HiFi训练的模型生成ONT数据 (及反向), 与同平台基线对比
        组合: ont→ont (基线), hifi→hifi (基线), hifi→ont (跨平台), ont→hifi (跨平台)
        """
        print("=" * 60)
        print("[消融实验3] 跨平台迁移验证")
        print("=" * 60)
        
        ref_fasta = self.data_dir / f'references/{test_species}_ref.fa'
        if not ref_fasta.exists():
            print(f"  ⚠ 参考基因组不存在: {ref_fasta}")
            return None
        
        combos = [('ont', 'ont'), ('hifi', 'hifi'), ('hifi', 'ont'), ('ont', 'hifi')]
        results = {}
        
        for train_plat, test_plat in combos:
            tag = f'{train_plat}_to_{test_plat}'
            bam_file = self.data_dir / f'real_reads_verified/{test_species}_{train_plat}_aligned.bam'
            real_fastq = self.data_dir / f'real_reads_verified/{test_species}_{test_plat}.fastq.gz'
            
            if not bam_file.exists() or not real_fastq.exists():
                print(f"\n--- {tag}: 数据缺失, 跳过 ---")
                results[tag] = {'note': 'real data not available'}
                continue
            
            print(f"\n--- {tag} ---")
            
            # 训练 (ONT基线模型通常已存在: {species}_errhmm.json)
            if train_plat == 'ont' and self._is_current_model(
                self.data_dir / f'trained_models/{test_species}_errhmm.json'
            ):
                model_path = self.data_dir / f'trained_models/{test_species}_errhmm.json'
            else:
                model_path = self.data_dir / f'trained_models/{test_species}_errhmm_{train_plat}{self.run_suffix}.json'
                if not self._is_current_model(model_path):
                    hmm = GCAwareErrHMM(gc_bins=10, window_size=100)
                    hmm.train_from_bam(
                        str(bam_file),
                        str(ref_fasta),
                        min_mapq=20,
                        max_reads=self.train_max_reads,
                    )
                    hmm.save(str(model_path))
            
            # 按目标平台参数生成模拟数据
            sim_prefix = self.data_dir / f'simulated/ablation/{test_species}_cross_{tag}{self.run_suffix}'
            sim_prefix.parent.mkdir(parents=True, exist_ok=True)
            
            simulator = GCAwareSimulator(
                errhmm_model=str(model_path),
                platform=test_plat,
                seed=self.seed,
            )
            sim_fastq = simulator.simulate_reads(
                str(ref_fasta), coverage=10,
                output_prefix=str(sim_prefix), threads=self.threads
            )
            
            # 与目标平台的真实数据对比
            eval_result = compare_real_vs_simulated(
                str(real_fastq), sim_fastq, str(ref_fasta)
            )
            results[tag] = {
                'composite_score': eval_result.get('composite_score', {}).get('composite_score', 0),
                'sub_scores': eval_result.get('composite_score', {}).get('sub_scores', {})
            }
            print(f"  {tag}: composite = {results[tag]['composite_score']}")
        
        # 迁移损失 = 同平台基线 - 跨平台得分
        for cross, base in [('hifi_to_ont', 'ont_to_ont'), ('ont_to_hifi', 'hifi_to_hifi')]:
            if cross in results and base in results and 'composite_score' in results[cross] and 'composite_score' in results[base]:
                results[f'{cross}_transfer_loss'] = round(
                    results[base]['composite_score'] - results[cross]['composite_score'], 2)
        
        if output_json:
            with open(output_json, 'w') as f:
                json.dump(results, f, indent=2)
        
        print(f"\n✓ 跨平台消融完成")
        return results
    
    def ablation_state_space(self, test_species='Ecoli', output_json=None):
        """
        消融实验4: 状态空间设计
        对比简化版(M/S/I/D) vs 完整版(M/S/I/D1/D2/D3/D4+)
        """
        print("=" * 60)
        print("[消融实验4] 状态空间设计对比")
        print("=" * 60)
        
        # 简化版：合并所有缺失长度为单一D状态
        # 完整版：区分D1, D2, D3, D4+
        
        ref_fasta = self.data_dir / f'references/{test_species}_ref.fa'
        bam_file = self.data_dir / f'real_reads_verified/{test_species}_ont_aligned.bam'
        
        if not bam_file.exists():
            print("  ⚠ BAM文件不存在")
            return None
        
        results = {}
        
        # 完整版 (默认)
        print("\n--- Full state space (M/S/I/D1/D2/D3/D4+) ---")
        hmm_full = GCAwareErrHMM(gc_bins=10, window_size=100)
        hmm_full.STATES = ['M', 'S', 'I', 'D1', 'D2', 'D3', 'D4+']
        hmm_full.train_from_bam(
            str(bam_file),
            str(ref_fasta),
            min_mapq=20,
            max_reads=self.train_max_reads,
        )
        
        model_full = self.data_dir / f'trained_models/{test_species}_errhmm_full{self.run_suffix}.json'
        hmm_full.save(str(model_full))
        
        # 简化版
        print("\n--- Simplified state space (M/S/I/D) ---")
        hmm_simple = GCAwareErrHMM(gc_bins=10, window_size=100)
        hmm_simple.STATES = ['M', 'S', 'I', 'D']
        hmm_simple.train_from_bam(
            str(bam_file),
            str(ref_fasta),
            min_mapq=20,
            max_reads=self.train_max_reads,
        )
        
        model_simple = self.data_dir / f'trained_models/{test_species}_errhmm_simple{self.run_suffix}.json'
        hmm_simple.save(str(model_simple))
        
        real_fastq = self.data_dir / f'real_reads_verified/{test_species}_ont.fastq.gz'
        profile_json = self.data_dir / f'tmp/read_profiles/{test_species}.json'

        def score_state_model(tag, model_path):
            sim_prefix = (
                self.data_dir / f'simulated/ablation/{test_species}_state_{tag}{self.run_suffix}'
            )
            simulator = GCAwareSimulator(
                errhmm_model=str(model_path),
                platform='ont',
                profile_json=str(profile_json) if profile_json.exists() else None,
                seed=self.seed,
            )
            sim_fastq = simulator.simulate_reads(
                str(ref_fasta), coverage=10,
                output_prefix=str(sim_prefix), threads=self.threads,
            )
            evaluation = compare_real_vs_simulated(
                str(real_fastq), sim_fastq, str(ref_fasta),
            )
            composite = evaluation.get('composite_score', {})
            return {
                'num_states': 7 if tag == 'full' else 4,
                'states': (
                    ['M', 'S', 'I', 'D1', 'D2', 'D3', 'D4+']
                    if tag == 'full' else ['M', 'S', 'I', 'D']
                ),
                'composite_score': composite.get('composite_score', 0),
                'sub_scores': composite.get('sub_scores', {}),
            }

        results['full_states'] = score_state_model('full', model_full)
        results['simple_states'] = score_state_model('simple', model_simple)
        if results['full_states']['composite_score'] >= results['simple_states']['composite_score']:
            results['best_state_space'] = 'full'
        else:
            results['best_state_space'] = 'simple'
        results['conclusion'] = (
            f"{results['best_state_space']} state space provides the higher "
            "Level-1 composite score"
        )
        
        if output_json:
            with open(output_json, 'w') as f:
                json.dump(results, f, indent=2)
        
        print(f"\n✓ 状态空间消融完成")
        return results


def main():
    parser = argparse.ArgumentParser(description='Ablation Study')
    parser.add_argument('study', choices=['gc-bins', 'coverage', 'cross-platform', 'state-space', 'all'])
    parser.add_argument('--project-dir', default='.')
    parser.add_argument('--species', default='Ecoli')
    parser.add_argument('--bins', type=int, default=10)
    parser.add_argument('--coverage', type=int)
    parser.add_argument('--output')
    parser.add_argument('--threads', type=int, default=16)
    
    args = parser.parse_args()
    
    study = AblationStudy(args.project_dir, threads=args.threads)
    
    if args.study == 'gc-bins':
        study.ablation_gc_bins(args.species, args.output)
    elif args.study == 'coverage':
        study.ablation_coverage(args.species, args.output)
    elif args.study == 'cross-platform':
        study.ablation_cross_platform(output_json=args.output)
    elif args.study == 'state-space':
        study.ablation_state_space(args.species, args.output)
    elif args.study == 'all':
        study.ablation_gc_bins(args.species, args.output.replace('.json', '_gc_bins.json') if args.output else None)
        study.ablation_coverage(args.species, args.output.replace('.json', '_coverage.json') if args.output else None)
        study.ablation_cross_platform(output_json=args.output.replace('.json', '_cross_platform.json') if args.output else None)
        study.ablation_state_space(args.species, args.output.replace('.json', '_state_space.json') if args.output else None)


if __name__ == '__main__':
    main()
