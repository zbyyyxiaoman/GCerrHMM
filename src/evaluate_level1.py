#!/usr/bin/env python3
"""
Level-1 质量评估模块
功能: 对比模拟数据与真实数据的基础质量指标

评估维度:
1. Read长度分布 (N50, mean, std)
2. 错误率 (total, substitution, insertion, deletion)
3. QV分布相关性 (Pearson)
4. GC含量分布 (Kolmogorov-Smirnov test)
5. k-mer频谱 (Jaccard similarity)
"""

import sys
import gzip
import json
import os
import random
import shutil
import subprocess
import tempfile
import numpy as np
import pysam
from Bio import SeqIO
from scipy import stats
from collections import Counter
import argparse
from pathlib import Path
import pickle

try:
    sys.stdout.reconfigure(line_buffering=True)
except Exception:
    pass


def _open_fastq(path):
    """按扩展名以文本模式打开 fastq (支持 .gz)"""
    if str(path).endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'rt')


class Level1Evaluator:
    """Level-1 模拟数据质量评估器"""
    
    def __init__(self, kmer_size=21, min_mapq=20):
        self.kmer_size = kmer_size
        self.max_reads = 5000
        # A full 21-mer Counter over thousands of long reads can exceed the
        # 7GB WSL memory limit. 250 reads still provides millions of k-mers
        # while keeping the comparison bounded.
        self.kmer_max_reads = 250
        self.kmer_sample_size = 10000
        self.min_mapq = min_mapq

    def analyze_fastq(self, fastq_path):
        """Single-pass sampling stats for long-read FASTQ files."""
        lengths = []
        qv_hist = np.zeros(94, dtype=np.int64)
        gc_values = []
        kmers = Counter()
        total_kmers = 0
        kmer_reads = 0
        window_size = 100

        rng = random.Random(20260913)
        reservoir = []
        for index, record in enumerate(SeqIO.parse(_open_fastq(fastq_path), 'fastq')):
            if index < self.max_reads:
                reservoir.append(record)
                continue
            slot = rng.randrange(index + 1)
            if slot < self.max_reads:
                reservoir[slot] = record

        for record in reservoir:
            seq = str(record.seq).upper()
            seq_len = len(seq)
            lengths.append(seq_len)

            qvs = record.letter_annotations.get('phred_quality')
            if qvs:
                qv_arr = np.asarray(qvs, dtype=np.int64)
                n_qv = min(94, qv_arr.size)
                qv_hist[:n_qv] += np.bincount(qv_arr, minlength=94)[:n_qv]

            for j in range(0, seq_len - window_size + 1, window_size):
                window = seq[j:j + window_size]
                gc_values.append((window.count('G') + window.count('C')) / window_size * 100)

            if kmer_reads < self.kmer_max_reads and seq_len >= self.kmer_size:
                for j in range(seq_len - self.kmer_size + 1):
                    kmer = seq[j:j + self.kmer_size]
                    if 'N' not in kmer:
                        kmers[kmer] += 1
                        total_kmers += 1
                kmer_reads += 1

        if not lengths:
            raise RuntimeError(f"No reads parsed from {fastq_path}")

        lengths_arr = np.asarray(lengths, dtype=float)
        read_hist = np.histogram(lengths_arr, bins=50)

        gc_arr = np.asarray(gc_values, dtype=float) if gc_values else np.zeros(1)
        gc_hist = np.histogram(gc_arr, bins=50, range=(0, 100))

        qv_total = int(qv_hist.sum())
        if qv_total > 0:
            idx = np.arange(94, dtype=float)
            cum = np.cumsum(qv_hist)
            qv_mean = float(np.average(idx, weights=qv_hist))
            qv_median = int(idx[np.searchsorted(cum, qv_total / 2)])
            qv_std = float(np.sqrt(np.average((idx - qv_mean) ** 2, weights=qv_hist)))
        else:
            qv_mean = qv_median = qv_std = 0.0

        qv_counts, _ = np.histogram(np.arange(94), weights=qv_hist, bins=50, range=(0, 93))

        return {
            'read_length': {
                'n50': self._compute_n50(lengths_arr),
                'mean': float(np.mean(lengths_arr)),
                'median': float(np.median(lengths_arr)),
                'std': float(np.std(lengths_arr)),
                'min': int(np.min(lengths_arr)),
                'max': int(np.max(lengths_arr)),
                'q25': float(np.percentile(lengths_arr, 25)),
                'q75': float(np.percentile(lengths_arr, 75)),
                'histogram': [h.tolist() for h in read_hist]
            },
            'qv': {
                'mean': qv_mean,
                'median': qv_median,
                'std': qv_std,
                'histogram': [qv_counts.tolist(), np.linspace(0, 93, 51).tolist()]
            },
            'gc': {
                'mean': float(np.mean(gc_arr)),
                'median': float(np.median(gc_arr)),
                'std': float(np.std(gc_arr)),
                'histogram': [h.tolist() for h in gc_hist]
            },
            'kmer': {
                'unique_kmers': len(kmers),
                'total_kmers': total_kmers,
                'kmer_freq': kmers
            },
            'n_reads_analyzed': len(lengths_arr),
            'n_kmer_reads_analyzed': kmer_reads
        }
    
    def _write_sampled_fastq(self, fastq_path, out_path, max_reads):
        """Write a deterministic reservoir sample of FASTQ records.

        Taking the first N records is not safe for SRA/ENA runs: project
        authors can prepend short, adapter-rich, or otherwise unusual reads,
        which biases the error-rate estimate.
        """
        rng = random.Random(20260913)
        reservoir = []
        for index, record in enumerate(SeqIO.parse(_open_fastq(fastq_path), 'fastq')):
            if index < max_reads:
                reservoir.append(record)
                continue
            slot = rng.randrange(index + 1)
            if slot < max_reads:
                reservoir[slot] = record

        with open(out_path, 'w') as out_handle:
            for record in reservoir:
                SeqIO.write([record], out_handle, 'fastq')
        return len(reservoir)

    def align_and_compute_error_rate(self, fastq_path, ref_fasta, tmp_dir, threads=4, max_reads=None):
        """Align a FASTQ to the reference and return error-rate statistics."""
        if max_reads is None:
            max_reads = self.max_reads
        try:
            bam_path = os.path.join(tmp_dir, "aligned.bam")
            sampled_fastq = os.path.join(tmp_dir, "sampled.fastq")
            sampled_reads_written = self._write_sampled_fastq(
                fastq_path,
                sampled_fastq,
                max_reads,
            )
            align_cmd = (
                f"minimap2 -ax map-ont -t {threads} -I 1G {ref_fasta} {sampled_fastq} | "
                f"samtools sort -@ {threads} -m 200M -o {bam_path}"
            )
            subprocess.run(
                align_cmd,
                shell=True,
                check=True,
                capture_output=True,
                text=True,
                timeout=1800,
            )
            subprocess.run(
                f"samtools index {bam_path}",
                shell=True,
                check=True,
                capture_output=True,
                text=True,
                timeout=600,
            )
            result = self.compute_error_rate(bam_path, ref_fasta)
            result['sampling_mode'] = 'deterministic_reservoir'
            result['sampled_reads_requested'] = int(max_reads)
            result['sampled_reads_written'] = int(sampled_reads_written)
            return result
        except Exception as exc:
            print(f"  [WARN] error-rate alignment failed: {exc}")
            return None
    
    # ========== 1. Read长度分布 ==========
    
    def compute_read_length_stats(self, fastq_path, max_reads=100000):
        """计算Read长度统计"""
        lengths = []
        for i, record in enumerate(SeqIO.parse(_open_fastq(fastq_path), 'fastq')):
            if i >= max_reads:
                break
            lengths.append(len(record.seq))
        
        lengths = np.array(lengths)
        
        return {
            'n50': self._compute_n50(lengths),
            'mean': float(np.mean(lengths)),
            'median': float(np.median(lengths)),
            'std': float(np.std(lengths)),
            'min': int(np.min(lengths)),
            'max': int(np.max(lengths)),
            'q25': float(np.percentile(lengths, 25)),
            'q75': float(np.percentile(lengths, 75)),
            'histogram': [h.tolist() for h in np.histogram(lengths, bins=50)]
        }
    
    def _compute_n50(self, lengths):
        """计算N50"""
        sorted_lengths = np.sort(lengths)[::-1]
        cumsum = np.cumsum(sorted_lengths)
        total = cumsum[-1]
        idx = np.searchsorted(cumsum, total / 2)
        return int(sorted_lengths[idx])
    
    # ========== 2. 错误率分析 ==========
    
    def compute_error_rate(self, bam_path, ref_fasta, max_reads=50000):
        """
        通过高质量 primary alignment 计算错误率。

        原始实现把 secondary/supplementary、MAPQ=0 以及几乎全部软裁剪的
        记录也计入分母，长读数据中少量嵌合/重复比对会把错误率放大到
        接近 100%。这里统一过滤 records，并仅统计有效比对区域。
        """
        ref_dict = {}
        for record in SeqIO.parse(ref_fasta, 'fasta'):
            ref_dict[record.id] = str(record.seq)
        
        bam = pysam.AlignmentFile(bam_path, 'rb')
        
        total_bases = 0
        mismatch = 0
        insertions = 0
        deletions = 0
        
        read_count = 0
        accepted_reads = 0
        mapq_sum = 0
        for read in bam:
            if (
                read.is_unmapped
                or read.is_secondary
                or read.is_supplementary
                or read.is_qcfail
                or read.is_duplicate
                or not read.cigartuples
            ):
                continue
            if read.mapping_quality < self.min_mapq:
                continue
            
            ref_name = read.reference_name
            if ref_name not in ref_dict:
                continue
            
            ref_seq = ref_dict[ref_name]
            # Pysam exposes the aligned query in reference orientation here.
            # Using query_sequence directly would reverse-strand reads appear
            # as ~75% mismatches and inflate the aggregate error rate.
            query_seq = read.query_alignment_sequence
            if query_seq is None:
                continue

            aligned_query_bases = sum(
                length
                for op, length in read.cigartuples
                if op in (0, 1, 7, 8)
            )
            if aligned_query_bases < 100:
                continue

            read_count += 1
            if read_count > max_reads:
                break
            accepted_reads += 1
            mapq_sum += read.mapping_quality
            ref_pos = read.reference_start
            query_pos = 0
            
            for cigar_op, length in read.cigartuples:
                if cigar_op in (0, 7, 8):  # M, =, X
                    for i in range(length):
                        if query_pos + i < len(query_seq):
                            ref_base = ref_seq[ref_pos + i].upper() if ref_pos + i < len(ref_seq) else 'N'
                            q_base = query_seq[query_pos + i].upper()
                            if ref_base != q_base:
                                mismatch += 1
                            total_bases += 1
                    ref_pos += length
                    query_pos += length
                    
                elif cigar_op == 1:  # I
                    insertions += length
                    total_bases += length
                    query_pos += length
                    
                elif cigar_op == 2:  # D
                    deletions += length
                    ref_pos += length
        
        bam.close()
        
        if total_bases == 0:
            return {
                'total': 0,
                'substitution': 0,
                'insertion': 0,
                'deletion': 0,
                'total_bases': 0,
                'accepted_reads': 0,
                'mean_mapq': 0,
                'min_mapq': self.min_mapq,
            }
        
        return {
            'total': (mismatch + insertions + deletions) / total_bases,
            'substitution': mismatch / total_bases,
            'insertion': insertions / total_bases,
            'deletion': deletions / total_bases,
            'total_bases': int(total_bases),
            'mismatches': int(mismatch),
            'inserted_bases': int(insertions),
            'deleted_bases': int(deletions),
            'accepted_reads': int(accepted_reads),
            'mean_mapq': mapq_sum / accepted_reads if accepted_reads else 0.0,
            'min_mapq': self.min_mapq,
        }
    
    # ========== 3. QV分布相关性 ==========
    
    def compute_qv_distribution(self, fastq_path, max_reads=100000):
        """计算质量值分布"""
        qvs = []
        for i, record in enumerate(SeqIO.parse(_open_fastq(fastq_path), 'fastq')):
            if i >= max_reads:
                break
            qvs.extend(record.letter_annotations['phred_quality'])
        
        qvs = np.array(qvs)
        
        return {
            'mean': float(np.mean(qvs)),
            'median': float(np.median(qvs)),
            'std': float(np.std(qvs)),
            'histogram': [h.tolist() for h in np.histogram(qvs, bins=50, range=(0, 93))]
        }
    
    def qv_correlation(self, real_qv, sim_qv):
        """计算QV分布的Pearson相关性"""
        # 使用直方图分布计算相关性
        real_hist = np.array(real_qv['histogram'][0], dtype=float)
        sim_hist = np.array(sim_qv['histogram'][0], dtype=float)
        
        # 归一化
        real_hist /= real_hist.sum() if real_hist.sum() > 0 else 1
        sim_hist /= sim_hist.sum() if sim_hist.sum() > 0 else 1
        
        if len(real_hist) != len(sim_hist):
            min_len = min(len(real_hist), len(sim_hist))
            real_hist = real_hist[:min_len]
            sim_hist = sim_hist[:min_len]
        
        corr, p_value = stats.pearsonr(real_hist, sim_hist)
        return {'pearson_r': float(corr), 'p_value': float(p_value)}
    
    # ========== 4. GC含量分布 ==========
    
    def compute_gc_distribution(self, fastq_path, window_size=100, max_reads=100000):
        """计算Read的GC含量分布"""
        gc_values = []
        for i, record in enumerate(SeqIO.parse(_open_fastq(fastq_path), 'fastq')):
            if i >= max_reads:
                break
            seq = str(record.seq).upper()
            for j in range(0, len(seq) - window_size + 1, window_size):
                window = seq[j:j + window_size]
                gc = (window.count('G') + window.count('C')) / len(window) * 100
                gc_values.append(gc)
        
        gc_values = np.array(gc_values)
        
        return {
            'mean': float(np.mean(gc_values)),
            'median': float(np.median(gc_values)),
            'std': float(np.std(gc_values)),
            'histogram': [h.tolist() for h in np.histogram(gc_values, bins=50, range=(0, 100))]
        }
    
    def gc_ks_test(self, real_gc, sim_gc):
        """GC分布的Kolmogorov-Smirnov检验"""
        real_samples = np.random.choice(
            np.linspace(0, 100, len(real_gc['histogram'][0])),
            size=10000,
            p=np.array(real_gc['histogram'][0]) / sum(real_gc['histogram'][0])
        )
        sim_samples = np.random.choice(
            np.linspace(0, 100, len(sim_gc['histogram'][0])),
            size=10000,
            p=np.array(sim_gc['histogram'][0]) / sum(sim_gc['histogram'][0])
        )
        
        ks_stat, p_value = stats.ks_2samp(real_samples, sim_samples)
        return {'ks_statistic': float(ks_stat), 'p_value': float(p_value)}
    
    # ========== 5. k-mer频谱 ==========
    
    def compute_kmer_spectrum(self, fastq_path, k=21, max_reads=50000):
        """计算k-mer频谱"""
        kmers = Counter()
        total_kmers = 0
        
        for i, record in enumerate(SeqIO.parse(_open_fastq(fastq_path), 'fastq')):
            if i >= max_reads:
                break
            seq = str(record.seq).upper()
            for j in range(len(seq) - k + 1):
                kmer = seq[j:j + k]
                if 'N' not in kmer:
                    kmers[kmer] += 1
                    total_kmers += 1
        
        # 转换为频率
        kmer_freq = {k: v / total_kmers for k, v in kmers.items()}
        
        return {
            'unique_kmers': len(kmers),
            'total_kmers': total_kmers,
            'kmer_freq': kmer_freq
        }
    
    def kmer_jaccard(self, real_kmer, sim_kmer):
        """计算k-mer频谱的Jaccard相似性"""
        real_freq = self._kmer_frequency_map(real_kmer)
        sim_freq = self._kmer_frequency_map(sim_kmer)
        real_set = set(real_freq.keys())
        sim_set = set(sim_freq.keys())
        
        intersection = len(real_set & sim_set)
        union = len(real_set | sim_set)
        
        jaccard = intersection / union if union > 0 else 0
        
        # 也计算频率加权相似度
        common_kmers = real_set & sim_set
        if len(common_kmers) > 0:
            real_freqs = np.array([real_freq[k] for k in common_kmers])
            sim_freqs = np.array([sim_freq[k] for k in common_kmers])
            freq_corr, _ = stats.pearsonr(real_freqs, sim_freqs)
        else:
            freq_corr = 0
        
        return {
            'jaccard': float(jaccard),
            'frequency_correlation': float(freq_corr)
        }

    def _kmer_frequency_map(self, kmer):
        """Return a normalized, comparably sampled frequency map."""
        if 'kmer_freq' in kmer:
            raw = kmer['kmer_freq']
            if hasattr(raw, 'most_common'):
                items = raw.most_common(self.kmer_sample_size)
                total = kmer.get('total_kmers', 0) or sum(raw.values()) or 1
                return {key: value / total for key, value in items}
            if isinstance(raw, dict):
                total = kmer.get('total_kmers', 0) or sum(raw.values()) or 1
                items = sorted(raw.items(), key=lambda item: item[1], reverse=True)
                return {
                    key: value / total
                    for key, value in items[:self.kmer_sample_size]
                }
        return kmer.get('kmer_freq_sample', {})
    
    # ========== 综合评分 ==========
    
    def compute_composite_score(self, metrics):
        """
        计算综合评分 (0-100)
        各维度权重可以根据需求调整
        """
        scores = {}
        if 'read_length_similarity' in metrics:
            scores['read_length'] = max(0.0, min(1.0, metrics['read_length_similarity'])) * 100
        if 'error_rate_similarity' in metrics:
            scores['error_rate'] = max(0.0, min(1.0, metrics['error_rate_similarity'])) * 100
        if 'qv_pearson' in metrics:
            scores['qv'] = max(0.0, min(1.0, (metrics['qv_pearson'] + 1.0) / 2.0)) * 100
        if 'gc_ks_stat' in metrics:
            scores['gc'] = max(0.0, 1.0 - metrics['gc_ks_stat']) * 100
        if metrics.get('kmer_freq_corr') is not None:
            scores['kmer'] = max(0.0, min(1.0, metrics['kmer_freq_corr'])) * 100
        elif 'kmer_jaccard' in metrics:
            scores['kmer'] = metrics['kmer_jaccard'] * 100

        weights = {
            'read_length': 0.15,
            'error_rate': 0.25,
            'qv': 0.20,
            'gc': 0.25,
            'kmer': 0.15
        }

        used_weights = sum(weights[k] for k in scores)
        composite = sum(scores[k] * weights[k] for k in scores)
        if used_weights > 0:
            composite = composite / used_weights

        return {
            'composite_score': round(composite, 2),
            'sub_scores': {k: round(v, 2) for k, v in scores.items()}
        }


def compare_real_vs_simulated(
    real_fastq,
    sim_fastq,
    ref_fasta=None,
    output_json=None,
    max_reads=None,
    kmer_max_reads=None,
    kmer_sample_size=None,
    cache_dir=None,
):
    """
    对比真实数据和模拟数据的Level-1指标
    
    Args:
        real_fastq: 真实数据FASTQ路径
        sim_fastq: 模拟数据FASTQ路径
        ref_fasta: 参考基因组 (用于错误率计算，可选)
        output_json: 输出JSON路径
    
    Returns:
        dict: 对比结果
    """
    print(f"[Level-1评估] 开始对比...")
    print(f"  真实数据: {real_fastq}")
    print(f"  模拟数据: {sim_fastq}")
    
    evaluator = Level1Evaluator()
    if max_reads is not None:
        evaluator.max_reads = int(max_reads)
    if kmer_max_reads is not None:
        evaluator.kmer_max_reads = int(kmer_max_reads)
    if kmer_sample_size is not None:
        evaluator.kmer_sample_size = int(kmer_sample_size)
    results = {'real': {}, 'simulated': {}, 'comparison': {}}

    # 1. Read长度分布
    print("  [1/5] 抽样统计真实/模拟数据...")
    cache_version = 5
    sampling_parameters = {
        'max_reads': evaluator.max_reads,
        'kmer_max_reads': evaluator.kmer_max_reads,
        'kmer_sample_size': evaluator.kmer_sample_size,
    }
    cache_dir = (
        Path(cache_dir).resolve()
        if cache_dir
        else Path(tempfile.gettempdir()) / 'gcerrhmm_level1_cache'
    )
    cache_dir.mkdir(parents=True, exist_ok=True)
    real_cache = cache_dir / (Path(real_fastq).name + '.json')
    real_stat = os.stat(real_fastq)
    if real_cache.exists():
        with open(real_cache) as f:
            cached = json.load(f)
        metadata = cached.get('_cache_metadata', {})
        if (
            metadata.get('version') == cache_version
            and metadata.get('source_size') == real_stat.st_size
            and metadata.get('source_mtime_ns') == real_stat.st_mtime_ns
            and metadata.get('max_reads') == sampling_parameters['max_reads']
            and metadata.get('kmer_max_reads') == sampling_parameters['kmer_max_reads']
            and metadata.get('kmer_sample_size') == sampling_parameters['kmer_sample_size']
        ):
            results['real'] = cached['stats']
        else:
            real_cache = None
    else:
        real_cache = None

    if real_cache is None:
        results['real'] = evaluator.analyze_fastq(real_fastq)
        kmer = results['real'].get('kmer', {})
        if 'kmer_freq' in kmer:
            counts = kmer['kmer_freq']
            if hasattr(counts, 'most_common'):
                top = counts.most_common(evaluator.kmer_sample_size)
            else:
                top = sorted(
                    counts.items(),
                    key=lambda kv: kv[1],
                    reverse=True,
                )[:evaluator.kmer_sample_size]
            total = kmer.get('total_kmers', 0) or 1
            kmer['kmer_freq_sample'] = {k: round(v / total, 8) for k, v in top}
            del kmer['kmer_freq']
        cache_path = cache_dir / (Path(real_fastq).name + '.json')
        temporary_cache = cache_path.with_name(
            cache_path.name + f'.tmp.{os.getpid()}'
        )
        with open(temporary_cache, 'w') as f:
            json.dump({
                '_cache_metadata': {
                    'version': cache_version,
                    'source_path': str(real_fastq),
                    'source_size': real_stat.st_size,
                    'source_mtime_ns': real_stat.st_mtime_ns,
                    'sampling': 'deterministic_reservoir',
                    **sampling_parameters,
                },
                'stats': results['real'],
            }, f, indent=2)
        os.replace(temporary_cache, cache_path)
    results['simulated'] = evaluator.analyze_fastq(sim_fastq)
    real_n50 = results['real']['read_length']['n50']
    sim_n50 = results['simulated']['read_length']['n50']
    results['comparison']['read_length_similarity'] = (
        min(real_n50, sim_n50) / max(real_n50, sim_n50) if max(real_n50, sim_n50) > 0 else 0.0
    )

    # 2. QV分布
    print("  [2/5] QV分布...")
    qv_corr = evaluator.qv_correlation(results['real']['qv'], results['simulated']['qv'])
    results['comparison']['qv_pearson'] = qv_corr['pearson_r']

    # 3. GC含量分布
    print("  [3/5] GC含量分布...")
    gc_ks = evaluator.gc_ks_test(results['real']['gc'], results['simulated']['gc'])
    results['comparison']['gc_ks_stat'] = gc_ks['ks_statistic']

    # 4. k-mer频谱
    print("  [4/5] k-mer频谱...")
    kmer_jac = evaluator.kmer_jaccard(results['real']['kmer'], results['simulated']['kmer'])
    results['comparison']['kmer_jaccard'] = kmer_jac['jaccard']
    results['comparison']['kmer_freq_corr'] = kmer_jac['frequency_correlation']

    # 5. 错误率 (需要参考基因组)
    if ref_fasta and os.environ.get("ENABLE_LEVEL1_ERROR_RATE") == "1":
        print("  [5/5] 错误率...")
        with tempfile.TemporaryDirectory(prefix='level1_err_') as tmp_dir:
            real_err = evaluator.align_and_compute_error_rate(real_fastq, ref_fasta, tmp_dir)
            sim_err = evaluator.align_and_compute_error_rate(sim_fastq, ref_fasta, tmp_dir)
        if real_err is not None and sim_err is not None:
            denom = max(real_err['total'], sim_err['total'])
            results['comparison']['error_rate_similarity'] = (
                max(0.0, 1.0 - abs(real_err['total'] - sim_err['total']) / denom)
                if denom > 0 else 0.0
            )
            results['comparison']['real_error_rate'] = real_err['total']
            results['comparison']['sim_error_rate'] = sim_err['total']

    # 综合评分
    print("  [计算综合评分...]")
    composite = evaluator.compute_composite_score(results['comparison'])
    results['composite_score'] = composite

    # 保留少量k-mer样例，避免把数百万条k-mer写进JSON
    for side in ('real', 'simulated'):
        kmer = results[side].get('kmer', {})
        if 'kmer_freq' in kmer:
            counts = kmer['kmer_freq']
            if hasattr(counts, 'most_common'):
                top = counts.most_common(1000)
            else:
                top = sorted(counts.items(), key=lambda kv: kv[1], reverse=True)[:1000]
            total = kmer.get('total_kmers', 0) or 1
            kmer['kmer_freq_sample'] = {k: round(v / total, 8) for k, v in top}
            del kmer['kmer_freq']

    print(f"\n✓ Level-1评估完成!")
    print(f"  综合评分: {composite['composite_score']}/100")
    print(f"  子项评分:")
    for k, v in composite['sub_scores'].items():
        print(f"    {k}: {v:.2f}")

    if output_json:
        with open(output_json, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\n✓ 结果已保存: {output_json}")

    return results


def main():
    parser = argparse.ArgumentParser(description='Level-1 Quality Evaluation')
    parser.add_argument('--real', required=True, help='Real FASTQ')
    parser.add_argument('--sim', required=True, help='Simulated FASTQ')
    parser.add_argument('--ref', help='Reference FASTA (for error rate)')
    parser.add_argument('--output', help='Output JSON path')
    
    args = parser.parse_args()
    compare_real_vs_simulated(args.real, args.sim, args.ref, args.output)


if __name__ == '__main__':
    main()
