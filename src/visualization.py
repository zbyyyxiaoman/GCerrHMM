#!/usr/bin/env python3
"""
结果可视化与统计分析模块
功能: 生成论文所需的图表和统计摘要

输出:
  - 箱线图 (3条路线 × 6物种)
  - 热力图 (GC含量 × 错误率)
  - 雷达图 (多维度综合评分)
  - 折线图 (消融实验结果)
  - 综合评分卡表格
"""

import sys
import json
import gzip
import random
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # 无GUI环境
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
from scipy import stats
from Bio import SeqIO
from evaluate_level1 import _open_fastq

# 设置中文字体和样式
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['figure.figsize'] = (10, 6)
sns.set_style('whitegrid')

# BMC Bioinformatics 配色方案
COLORS = {
    'route_A_sample': '#3498db',  # 蓝
    'route_B_qshmm': '#e74c3c',   # 红
    'route_C_errhmm': '#2ecc71',  # 绿 (我们的方法)
    'highlight': '#f39c12',        # 橙
}

ROUTES = ['route_A_sample', 'route_B_qshmm', 'route_C_errhmm']
ROUTE_LABELS = ['Sample', 'qsHMM(nongc)', 'errHMM (Ours)']
# 与 config/species_config_v2.json 一致的六级物种梯度
SPECIES_ORDER = ['Ecoli', 'Scerevisiae', 'Athaliana', 'Dmelanogaster', 'Mmusculus_chr19', 'Hsapiens_chr21']
SPECIES_LABELS = ['E. coli\n(50.8% GC)', 'S. cerevisiae\n(38.2% GC)', 'A. thaliana\n(36.0% GC)',
                  'D. melanogaster\n(42.0% GC)', 'M. musculus\nchr19 (42.0% GC)', 'H. sapiens\nchr21 (40.9% GC)']
# 外部基线 (PBSIM3), 若存在也一并画
EXT_ROUTES = ['pbsim3_sample', 'pbsim3_qshmm']


def parse_result_key(key):
    """
    从结果文件key解析 (species, route)
    key形如 {species}_{route}_r{rep}, 例如 Ecoli_route_A_sample_r1
    """
    for route in ROUTES + EXT_ROUTES:
        if f'_{route}_' in key:
            return key.split(f'_{route}_')[0], route
    return None, None


class ResultVisualizer:
    """结果可视化器"""
    
    def __init__(self, results_dir, output_dir, models_dir='data/trained_models', data_dir='data'):
        self.results_dir = Path(results_dir)
        self.output_dir = Path(output_dir)
        self.models_dir = Path(models_dir)
        self.data_dir = Path(data_dir)
        self.real_dir = self.data_dir / 'real_reads_verified'
        self.sim_dir = self.data_dir / 'simulated'
        self.ref_dir = self.data_dir / 'reference'
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 加载所有结果
        self.level1_results = self._load_level1_results()
        self.mapping_results = self._load_mapping_results()
        self.variant_results = self._load_variant_results()
        self.assembly_results = self._load_assembly_results()
    
    def _load_level1_results(self):
        """加载Level-1评估结果"""
        results = {}
        for f in self.results_dir.glob('level1_*.json'):
            with open(f) as fp:
                key = f.stem.replace('level1_', '')
                results[key] = json.load(fp)
        return results
    
    def _load_mapping_results(self):
        """加载比对评估结果"""
        results = {}
        for f in self.results_dir.glob('mapping_*.json'):
            with open(f) as fp:
                key = f.stem.replace('mapping_', '')
                results[key] = json.load(fp)
        return results
    
    def _load_variant_results(self):
        """加载变异检测结果"""
        results = {}
        for f in self.results_dir.glob('variant_*.json'):
            with open(f) as fp:
                key = f.stem.replace('variant_', '')
                results[key] = json.load(fp)
        return results
    
    def _load_assembly_results(self):
        """加载拼接评估结果"""
        results = {}
        for f in self.results_dir.glob('assembly_*.json'):
            with open(f) as fp:
                key = f.stem.replace('assembly_', '')
                results[key] = json.load(fp)
        return results
    
    def plot_composite_score_boxplot(self):
        """综合评分箱线图 (Figure 1)"""
        print("[可视化] 生成综合评分箱线图...")
        
        # 整理数据
        data = []
        for key, result in self.level1_results.items():
            species, route = parse_result_key(key)
            if route is None:
                continue
            score = result.get('composite_score', {}).get('composite_score', 0)
            data.append({'Species': species, 'Route': route, 'Score': score})
        
        df = pd.DataFrame(data)
        if df.empty:
            print("  ⚠ 无Level-1数据, 跳过fig1")
            return
        
        # 只画实际有数据的物种
        species_present = [s for s in SPECIES_ORDER if s in set(df['Species'])]
        extra = sorted(set(df['Species']) - set(SPECIES_ORDER))
        species_present += extra
        
        # 创建图
        fig, ax = plt.subplots(figsize=(12, 7))
        
        positions = []
        box_data = []
        colors = []
        
        for i, species in enumerate(species_present):
            for j, route in enumerate(ROUTES):
                pos = i * 4 + j
                positions.append(pos)
                
                subset = df[(df['Species'] == species) & (df['Route'] == route)]['Score']
                box_data.append(subset.values if len(subset) > 0 else [0])
                colors.append(COLORS[route])
        
        bp = ax.boxplot(box_data, positions=positions, widths=0.6, 
                        patch_artist=True, showfliers=False)
        
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        # 设置x轴
        label_map = dict(zip(SPECIES_ORDER, SPECIES_LABELS))
        ax.set_xticks([i * 4 + 1 for i in range(len(species_present))])
        ax.set_xticklabels([label_map.get(s, s) for s in species_present], fontsize=9)
        ax.set_ylabel('Composite Quality Score', fontsize=12)
        ax.set_title('Level-1 Composite Score: 3 Methods x 6 Species (n=2 per cell)',
                     fontsize=14, fontweight='bold')
        ax.set_ylim(0, 105)
        
        # 图例
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=COLORS[r], alpha=0.7, label=l) 
                          for r, l in zip(ROUTES, ROUTE_LABELS)]
        ax.legend(handles=legend_elements, loc='lower left', fontsize=10)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'fig1_composite_score_boxplot.png')
        plt.close()
        print(f"  ✓ 已保存: fig1_composite_score_boxplot.png")
    
    def plot_gc_error_heatmap(self):
        """GC含量 vs 错误率热力图 (Figure 2)
        
        数据来源: 训练好的模型文件 (data/trained_models/), 非编造数据
          - Route C: errHMM各GC_bin的错误概率 1 - P(M→M)
          - 每个物种保留独立一行，不跨物种平均
        """
        print("[可视化] 生成GC-错误率热力图...")
        
        # Route C: 从errHMM模型提取 per-bin 错误率
        errhmm_curves = {}
        for model_file in sorted(self.models_dir.glob('*_errhmm.json')):
            species = model_file.name.removesuffix('_errhmm.json')
            with open(model_file) as f:
                model = json.load(f)
            states = model['states']
            m_idx = states.index('M')
            tp = model['transition_probs']
            per_bin = []
            for b in range(10):
                if str(b) in tp:
                    per_bin.append(1.0 - tp[str(b)][m_idx][m_idx])
            if len(per_bin) == 10:
                errhmm_curves[species] = per_bin
        
        if not errhmm_curves:
            print("  ⚠ 无errHMM模型数据, 跳过fig2")
            return

        species_order = sorted(errhmm_curves)
        data = np.vstack([
            np.asarray(errhmm_curves[species], dtype=float)
            / np.mean(errhmm_curves[species])
            for species in species_order
        ])
        
        fig, ax = plt.subplots(figsize=(10, max(4, 0.75 * len(species_order) + 1.5)))
        
        gc_labels = [f'{i*10}-{(i+1)*10}%' for i in range(10)]
        
        im = ax.imshow(data, cmap='RdYlGn_r', aspect='auto')
        
        ax.set_xticks(range(10))
        ax.set_xticklabels(gc_labels, rotation=45, ha='right')
        ax.set_yticks(range(len(species_order)))
        ax.set_yticklabels([species.replace('_', ' ') for species in species_order])
        
        ax.set_xlabel('GC Content Bin', fontsize=12)
        ax.set_title('errHMM Error Rate by GC Content (per species, relative to mean)',
                     fontsize=13, fontweight='bold')
        
        # 添加数值标注
        for i in range(len(species_order)):
            for j in range(10):
                ax.text(j, i, f'{data[i, j]:.2f}', ha='center', va='center', fontsize=8)
        
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Relative Error Rate', rotation=270, labelpad=20)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'fig2_gc_error_heatmap.png')
        plt.close()
        print(f"  ✓ 已保存: fig2_gc_error_heatmap.png (基于{len(species_order)}个物种模型)")
    
    def plot_quality_dimensions(self):
        """Figure 3: grouped bar chart of Level-1 quality dimensions."""
        print("[可视化] 生成质量维度条形图...")
        metric_keys = ['read_length', 'error_rate', 'qv', 'gc', 'kmer']
        category_labels = ['Read Length', 'Error Rate', 'QV Corr.', 'GC Match', 'k-mer Sim.']
        per_route = {route: {k: [] for k in metric_keys} for route in ROUTES}
        present_keys = set(metric_keys)

        for key, result in self.level1_results.items():
            _, route = parse_result_key(key)
            if route not in per_route:
                continue
            sub_scores = result.get('composite_score', {}).get('sub_scores', {})
            present_keys.intersection_update(
                k for k in metric_keys if sub_scores.get(k) is not None
            )
            for k in metric_keys:
                value = sub_scores.get(k)
                if value is not None:
                    per_route[route][k].append(value)

        keys = [k for k in metric_keys if k in present_keys]
        if not keys or not any(any(per_route[r][k]) for r in ROUTES for k in keys):
            print("  ⚠ 无可用Level-1维度, 跳过fig3")
            return

        x = np.arange(len(keys))
        width = 0.25
        fig, ax = plt.subplots(figsize=(11, 6))
        for i, (route, label) in enumerate(zip(ROUTES, ROUTE_LABELS)):
            means = []
            stds = []
            for k in keys:
                vals = per_route[route][k]
                means.append(float(np.mean(vals)) if vals else 0.0)
                stds.append(float(np.std(vals)) if len(vals) > 1 else 0.0)
            ax.bar(x + (i - 1) * width, means, width, yerr=stds, capsize=3,
                   label=label, color=COLORS[route], alpha=0.8)

        ax.set_xticks(x)
        ax.set_xticklabels(
            [label for k, label in zip(metric_keys, category_labels) if k in keys],
            fontsize=10,
        )
        ax.set_ylim(0, 100)
        ax.set_ylabel('Mean Score')
        ax.set_title('Level-1 Quality Dimensions by Route (mean ± SD)',
                     fontsize=13, fontweight='bold')
        ax.legend(fontsize=9, loc='lower left')
        ax.grid(axis='y', alpha=0.3)
        plt.tight_layout()
        plt.savefig(self.output_dir / 'fig3_quality_dimensions.png')
        plt.close()
        print(f"  ✓ 已保存: fig3_quality_dimensions.png")
    
    def plot_ablation_line(self):
        """消融实验折线图 (Figure 4)
        
        数据来源: results/tables/ablation_gcbins_*.json, ablation_cov_*.json
        无数据时跳过, 不使用任何硬编码数值
        """
        print("[可视化] 生成消融实验图...")
        
        # GC-bin粒度消融: 汇总所有 ablation_gcbins_*.json
        gc_series = {}
        stats_dir = self.results_dir.parent / 'stats'
        summary_files = sorted(
            stats_dir.glob('gc_bins_cross_species_with_1bin_*seeds*.csv')
        )
        if not summary_files:
            summary_files = sorted(
                stats_dir.glob('profile_matched_gc_bins_seeds_*.csv')
            )
        if summary_files:
            import csv
            with summary_files[-1].open() as handle:
                for row in csv.DictReader(handle):
                    species = row['species']
                    bins = int(row['gc_bins'].replace('bin', ''))
                    gc_series.setdefault(species, {})[bins] = (
                        float(row.get('mean', row.get('composite_score', 0))),
                        float(row.get('std', 0)),
                    )
        
        # 覆盖度消融: 优先使用最新的 3-repeat summary，避免混入旧 run 或
        # 将重复实验的最后一个文件误当成全量结果。
        cov_scores_dict = {}
        cov_sd_dict = {}
        coverage_summary_files = sorted(
            (self.results_dir.parent / 'stats').glob(
                'coverage_replicates_[0-9]*.csv'
            )
        )
        if coverage_summary_files:
            import csv
            with coverage_summary_files[-1].open() as handle:
                for row in csv.DictReader(handle):
                    coverage = int(row['coverage'].replace('x', ''))
                    cov_scores_dict[coverage] = float(row['mean'])
                    cov_sd_dict[coverage] = float(row['std'])
        else:
            for pattern in ('ablation_cov_*.json', 'ablation_cov.json'):
                for f in sorted(self.results_dir.glob(pattern)):
                    with open(f) as fp:
                        d = json.load(fp)
                    for k, v in d.items():
                        if k.endswith('x') and isinstance(v, dict) and 'composite_score' in v:
                            cov_scores_dict[int(k[:-1])] = v['composite_score']
        
        if not gc_series and not cov_scores_dict:
            print("  ⚠ 无消融实验数据, 跳过fig4")
            return
        
        coverages = sorted(cov_scores_dict)
        cov_scores = [cov_scores_dict[c] for c in coverages]
        cov_errors = [cov_sd_dict.get(c, 0.0) for c in coverages]
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
        
        # GC-bin消融
        if gc_series:
            all_bins = sorted({b for scores in gc_series.values() for b in scores})
            all_scores = [
                score[0] for scores in gc_series.values() for score in scores.values()
            ]
            palette = {
                'Ecoli': COLORS['route_A_sample'],
                'Athaliana': COLORS['route_C_errhmm'],
            }
            for species, scores in sorted(gc_series.items()):
                bins = sorted(scores)
                values = [scores[b][0] for b in bins]
                errors = [scores[b][1] for b in bins]
                ax1.errorbar(
                    bins, values, yerr=errors, fmt='o-', capsize=4,
                    color=palette.get(species, COLORS['route_B_qshmm']),
                    linewidth=2, markersize=8, label=species,
                )
                for b, s in zip(bins, values):
                    ax1.text(b, s + 0.25, f'{s:.1f}', ha='center', fontsize=8)
            ax1.set_xticks(all_bins)
            ax1.set_ylim(min(all_scores) - 2, max(all_scores) + 2)
            ax1.legend(fontsize=8)
        else:
            ax1.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax1.transAxes)
        ax1.set_xlabel('Number of GC Bins', fontsize=12)
        ax1.set_ylabel('Composite Score', fontsize=12)
        ax1.set_title('(A) GC-Bin Granularity', fontsize=13, fontweight='bold')
        ax1.grid(True, alpha=0.3)
        ax1.text(0.02, 0.02, 'Profile-matched, mean +/- SD across 3 seeds',
                 transform=ax1.transAxes, fontsize=8)
        
        # 覆盖度消融
        if coverages:
            ax2.errorbar(
                coverages,
                cov_scores,
                yerr=cov_errors,
                fmt='s-',
                capsize=4,
                color=COLORS['route_B_qshmm'],
                linewidth=2,
                markersize=8,
            )
            ax2.set_xticks(coverages)
            lower = min(score - error for score, error in zip(cov_scores, cov_errors))
            upper = max(score + error for score, error in zip(cov_scores, cov_errors))
            ax2.set_ylim(lower - 3, upper + 3)
            for c, s in zip(coverages, cov_scores):
                ax2.text(c, s + 0.5, f'{s:.1f}', ha='center', fontsize=9)
            ax2.text(
                0.02,
                0.02,
                'Mean +/- SD across 3 seeds',
                transform=ax2.transAxes,
                fontsize=8,
            )
        else:
            ax2.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax2.transAxes)
        ax2.set_xlabel('Training Coverage', fontsize=12)
        ax2.set_ylabel('Composite Score', fontsize=12)
        ax2.set_title('(B) Training Data Coverage', fontsize=13, fontweight='bold')
        ax2.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.output_dir / 'fig4_ablation_study.png')
        plt.close()
        print(f"  ✓ 已保存: fig4_ablation_study.png")
    
    def plot_downstream_comparison(self):
        """下游任务对比柱状图 (Figure 5)
        
        数据来源: results/tables/ 下真实的 mapping/variant/assembly/phasing 结果
        每个指标按路线取均值, 再按该指标的最大值归一化到0-100便于同图展示
        无数据时跳过, 不使用任何硬编码数值
        """
        print("[可视化] 生成下游任务对比图...")
        
        def mean_by_route(results_dict, metric_key, transform=lambda x: x):
            """按路线聚合某指标的均值"""
            per_route = {r: [] for r in ROUTES}
            for key, result in results_dict.items():
                _, route = parse_result_key(key)
                if route in per_route and metric_key in result:
                    per_route[route].append(transform(result[metric_key]))
            return {r: (float(np.mean(v)) if v else None) for r, v in per_route.items()}
        
        # 加载phasing结果 (此处按需加载, 避免__init__冗余)
        phasing_results = {}
        for f in self.results_dir.glob('phasing_*.json'):
            with open(f) as fp:
                phasing_results[f.stem.replace('phasing_', '')] = json.load(fp)
        
        metrics = {
            'Mapping\nRate':      mean_by_route(self.mapping_results, 'alignment_rate', lambda x: x * 100),
            'Variant\nF1':        mean_by_route(self.variant_results, 'f1_score', lambda x: x * 100),
            'Assembly\nN50':      mean_by_route(self.assembly_results, 'n50'),
            'Phasing\nAccuracy':  mean_by_route(phasing_results, 'switch_error_rate', lambda x: (1 - x) * 100),
        }
        # 只保留有数据的指标
        metrics = {k: v for k, v in metrics.items() if any(x is not None for x in v.values())}
        
        if not metrics:
            print("  ⚠ 无下游任务数据, 跳过fig5")
            return
        
        tasks = list(metrics.keys())
        # 每个指标按最大值归一化 (N50与其他指标量纲不同)
        series = {r: [] for r in ROUTES}
        for task in tasks:
            vals = metrics[task]
            vmax = max(v for v in vals.values() if v is not None) or 1
            for r in ROUTES:
                series[r].append(vals[r] / vmax * 100 if vals[r] is not None else 0)
        
        x = np.arange(len(tasks))
        width = 0.25
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        for i, (route, label) in enumerate(zip(ROUTES, ROUTE_LABELS)):
            ax.bar(x + (i - 1) * width, series[route], width, label=label,
                   color=COLORS[route], alpha=0.8)
        
        ax.set_ylabel('Relative Performance (% of best)', fontsize=12)
        ax.set_title('Downstream Task Performance Comparison', fontsize=14, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(tasks, fontsize=11)
        ax.set_ylim(0, 110)
        ax.legend(fontsize=9, loc='upper left', bbox_to_anchor=(1.01, 1))
        ax.grid(axis='y', alpha=0.3)
        
        plt.tight_layout(rect=[0, 0, 0.86, 1])
        plt.savefig(self.output_dir / 'fig5_downstream_comparison.png')
        plt.close()
        print(f"  ✓ 已保存: fig5_downstream_comparison.png")
    
    def generate_scorecard_table(self):
        """生成综合评分卡表格"""
        print("[可视化] 生成评分卡表格...")
        
        table_data = []
        
        # 按 (species, route) 聚合所有重复
        grouped = {}
        for key, result in self.level1_results.items():
            species, route = parse_result_key(key)
            if route is None:
                continue
            score = result.get('composite_score', {}).get('composite_score', 0)
            grouped.setdefault((species, route), []).append(score)
        
        label_map = dict(zip(ROUTES, ROUTE_LABELS))
        species_sorted = [s for s in SPECIES_ORDER 
                          if any(sp == s for sp, _ in grouped)]
        species_sorted += sorted({sp for sp, _ in grouped} - set(SPECIES_ORDER))
        
        for species in species_sorted:
            for route in ROUTES:
                scores = grouped.get((species, route))
                if scores:
                    table_data.append({
                        'Species': species,
                        'Method': label_map[route],
                        'Replicates': len(scores),
                        'Mean Score': f'{np.mean(scores):.2f}',
                        'Std': f'{np.std(scores):.2f}',
                        'Best': f'{np.max(scores):.2f}',
                        'Worst': f'{np.min(scores):.2f}'
                    })
        
        if not table_data:
            print("  ⚠ 无Level-1数据, 跳过评分卡")
            return
        
        df = pd.DataFrame(table_data)
        
        # 保存CSV
        csv_path = self.output_dir / 'scorecard_table.csv'
        df.to_csv(csv_path, index=False)
        
        # 保存Markdown表格
        md_path = self.output_dir / 'scorecard_table.md'
        with open(md_path, 'w') as f:
            f.write('# Composite Score Card\n\n')
            headers = list(df.columns)
            f.write('| ' + ' | '.join(headers) + ' |\n')
            f.write('|' + '|'.join(['---'] * len(headers)) + '|\n')
            for row in df.itertuples(index=False):
                f.write('| ' + ' | '.join(str(x) for x in row) + ' |\n')
        
        print(f"  ✓ 评分卡已保存: {csv_path}, {md_path}")
    
    def _species_real_fastq(self, species):
        path = self.real_dir / f"{species}_ont.fastq.gz"
        if path.exists():
            return path
        fallback = {
            'Hsapiens_chr21': 'Hsapiens_ont.fastq.gz',
            'Mmusculus_chr19': 'Mmusculus_ont.fastq.gz',
        }
        if species in fallback:
            path = self.real_dir / fallback[species]
            if path.exists():
                return path
        return None

    def _read_length_samples(self, path, max_reads=5000):
        rng = random.Random(20260913)
        reservoir = []
        for index, record in enumerate(SeqIO.parse(_open_fastq(path), 'fastq')):
            if index < max_reads:
                reservoir.append(len(record.seq))
                continue
            slot = rng.randrange(index + 1)
            if slot < max_reads:
                reservoir[slot] = len(record.seq)
        return reservoir

    def plot_pbsim3_read_length_distribution(self):
        """PBSIM3-style read length distribution by route and species."""
        print("[可视化] 生成 PBSIM3 风格 read length 分布图...")
        species = [s for s in SPECIES_ORDER if self._species_real_fastq(s) is not None]
        if not species:
            print("  [WARN] no real FASTQ found for read length figure")
            return

        fig, axes = plt.subplots(
            len(species), 1,
            figsize=(9, 3.2 * len(species)),
            sharex=True,
        )
        if len(species) == 1:
            axes = [axes]

        for ax, sp in zip(axes, species):
            real_path = self._species_real_fastq(sp)
            real_lengths = np.asarray(self._read_length_samples(real_path), dtype=float)
            ax.plot(
                np.sort(real_lengths),
                np.linspace(0, 1, len(real_lengths)),
                color='black',
                linewidth=2,
                label='Real ONT',
            )
            for route in ROUTES:
                vals = []
                for rep in (1, 2):
                    sim = self.sim_dir / route / sp / f"{sp}_{route}_r{rep}.fastq.gz"
                    if sim.exists():
                        vals.extend(self._read_length_samples(sim, max_reads=3000))
                if vals:
                    vals = np.sort(np.asarray(vals, dtype=float))
                    ax.plot(
                        vals,
                        np.linspace(0, 1, len(vals)),
                        color=COLORS[route],
                        linewidth=1.5,
                        label=ROUTE_LABELS[ROUTES.index(route)],
                    )
            ax.set_title(sp.replace('_', ' '))
            ax.set_xlim(0, 60000)
            ax.set_ylabel("CDF")
            ax.legend(fontsize=8, loc='lower right')

        axes[-1].set_xlabel("Read length (bp)")
        fig.suptitle("PBSIM3-style read length distribution: real vs simulated", y=1.02)
        fig.tight_layout()
        out = self.output_dir / 'fig6_pbsim3_read_length_distribution.png'
        fig.savefig(out, bbox_inches='tight')
        plt.close(fig)
        print(f"  ✓ 已保存: {out.name}")

    def _load_contig_lengths(self, species):
        ref = self.ref_dir / f"{species}_ref.fa"
        if not ref.exists():
            return {}
        return {r.id: len(r.seq) for r in SeqIO.parse(ref, 'fasta')}

    def plot_pbsim3_start_position_bias(self):
        """PBSIM3-style start-position uniformity summarized by KL divergence."""
        print("[可视化] 生成 PBSIM3 风格 read start position 统计图...")
        fig, axes = plt.subplots(
            len(SPECIES_ORDER), 1,
            figsize=(9, 3.0 * len(SPECIES_ORDER)),
            sharex=False,
        )
        used = False

        for ax, sp in zip(axes, SPECIES_ORDER):
            contig_lengths = self._load_contig_lengths(sp)
            total_ref = sum(contig_lengths.values())
            if total_ref <= 0:
                continue
            kl_by_route = {route: [] for route in ROUTES}
            for route in ROUTES:
                norm_vals = []
                for rep in (1, 2):
                    sim = self.sim_dir / route / sp / f"{sp}_{route}_r{rep}.fastq.gz"
                    if not sim.exists():
                        continue
                    reads_in_file = 0
                    for rec in SeqIO.parse(_open_fastq(sim), 'fastq'):
                        reads_in_file += 1
                        if reads_in_file > 5000:
                            break
                        desc = rec.description
                        pos_m = re.search(r'pos=(\d+)', desc)
                        contig_m = re.search(r'contig=(\S+)', desc)
                        if not pos_m:
                            continue
                        pos = int(pos_m.group(1))
                        contig = contig_m.group(1) if contig_m else None
                        length = contig_lengths.get(contig, total_ref)
                        if length <= 0:
                            continue
                        norm_vals.append(pos / length)
                if len(norm_vals) < 20:
                    continue
                hist, _ = np.histogram(norm_vals, bins=20, range=(0, 1))
                p = hist / hist.sum()
                q = np.full_like(p, 1.0 / len(p))
                eps = 1e-12
                kl_by_route[route].append(
                    float(np.sum(p * np.log2((p + eps) / q)))
                )

            if not any(kl_by_route.values()):
                continue
            x = np.arange(len(ROUTES))
            means = [
                float(np.mean(kl_by_route[r])) if kl_by_route[r] else 0.0
                for r in ROUTES
            ]
            stds = [
                float(np.std(kl_by_route[r])) if len(kl_by_route[r]) > 1 else 0.0
                for r in ROUTES
            ]
            ax.bar(x, means, yerr=stds, capsize=3,
                   color=[COLORS[r] for r in ROUTES], alpha=0.8)
            ax.set_xticks(x)
            ax.set_xticklabels(ROUTE_LABELS, fontsize=8)
            ax.set_title(sp.replace('_', ' '))
            ax.set_ylabel("KL vs uniform")
            ax.grid(axis='y', alpha=0.3)
            used = True

        if not used:
            plt.close(fig)
            print("  [WARN] no simulated start-position metadata found")
            return

        axes[-1].set_xlabel("Simulation route")
        fig.suptitle(
            "PBSIM3-style start-position uniformity (higher = less uniform)",
            y=1.02,
        )
        fig.tight_layout()
        out = self.output_dir / 'fig7_pbsim3_start_position_bias.png'
        fig.savefig(out, bbox_inches='tight')
        plt.close(fig)
        print(f"  ✓ 已保存: {out.name}")

    def generate_all(self):
        """生成所有图表"""
        print("\n===== 生成所有可视化 =====")
        self.plot_composite_score_boxplot()
        self.plot_gc_error_heatmap()
        self.plot_quality_dimensions()
        self.plot_ablation_line()
        self.plot_downstream_comparison()
        self.plot_pbsim3_read_length_distribution()
        self.plot_pbsim3_start_position_bias()
        self.generate_scorecard_table()
        print("\n✓ 所有可视化完成!")
        print(f"  输出目录: {self.output_dir}")


def main():
    parser = argparse.ArgumentParser(description='Generate Visualization and Statistics')
    parser.add_argument('--results-dir', required=True, help='Results JSON directory')
    parser.add_argument('--output-dir', required=True, help='Output figures directory')
    parser.add_argument('--stats-dir', help='Statistics output directory')
    parser.add_argument('--models-dir', default='data/trained_models',
                        help='Trained HMM models directory (for fig2)')
    parser.add_argument('--data-dir', default='data',
                        help='Project data directory')
    
    args = parser.parse_args()
    
    visualizer = ResultVisualizer(
        args.results_dir,
        args.output_dir,
        args.models_dir,
        args.data_dir,
    )
    visualizer.generate_all()


if __name__ == '__main__':
    main()
