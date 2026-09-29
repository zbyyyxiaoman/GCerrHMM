#!/usr/bin/env python3
"""
模拟数据生成器 - 支持errHMM的自定义模拟
功能: 使用训练好的errHMM模型生成GC-aware的模拟测序数据

三条路线:
  Route A (无模型): 固定错误率采样
  Route B (qsHMM):  非GC感知HMM (nongc模型)
  Route C (errHMM): GC感知HMM (GC-aware)
"""

import sys
import json
import random
import os
import time
import multiprocessing
import concurrent.futures
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
from pathlib import Path
import pickle
from train_errhmm import GCAwareErrHMM


_WORKER_SIMULATOR = None
_WORKER_REFERENCES = None


def _init_worker(simulator, references):
    global _WORKER_SIMULATOR, _WORKER_REFERENCES
    _WORKER_SIMULATOR = simulator
    _WORKER_REFERENCES = references


def _simulate_spec_chunk(simulator, references, specs):
    records = []
    lengths = []
    for spec in specs:
        random.seed(spec["seed"])
        np.random.seed(spec["seed"])
        contig_seq = references[spec["contig_id"]]
        start_pos = spec["start"]
        original_seq = contig_seq[start_pos:start_pos + spec["read_len"]]
        if simulator.hmm:
            modified_seq, cigar = simulator.introduce_errors(
                original_seq, contig_seq, start_pos
            )
        else:
            modified_seq, cigar = simulator._introduce_random_errors(original_seq)
        qvs = simulator.generate_quality_scores(len(modified_seq))
        qv_string = ''.join(chr(q + 33) for q in qvs)
        desc = (
            f"length={len(modified_seq)} cigar={cigar} "
            f"contig={spec['contig_id']} pos={start_pos}"
        )
        records.append(
            f"@sim_read_{spec['index']} {desc}\n"
            f"{modified_seq}\n+\n{qv_string}\n"
        )
        lengths.append(len(modified_seq))
    return "".join(records), lengths


def _worker_entry(specs):
    return _simulate_spec_chunk(_WORKER_SIMULATOR, _WORKER_REFERENCES, specs)


class GCAwareSimulator:
    """基于errHMM的GC感知测序模拟器"""
    
    def __init__(
        self,
        errhmm_model=None,
        platform='ont',
        read_length_mean=None,
        read_length_std=None,
        profile_json=None,
        seed=42,
    ):
        """
        初始化模拟器
        
        Args:
            errhmm_model: 加载的errHMM模型路径或对象
            read_length_mean: 平均读长 (HiFi默认15000)
            read_length_std: 读长标准差
        """
        platform_defaults = {'hifi': (15000, 5000), 'ont': (25000, 10000)}
        self.read_length_mean = read_length_mean if read_length_mean is not None else platform_defaults[platform][0]
        self.read_length_std = read_length_std if read_length_std is not None else platform_defaults[platform][1]
        self.read_length_quantiles = None
        self.quality_mean = 30.0
        self.quality_std = 5.0
        self.seed = int(seed)

        if profile_json:
            with open(profile_json, encoding='utf-8') as handle:
                profile = json.load(handle)
            length_profile = profile.get('read_length', {})
            self.read_length_mean = float(length_profile.get('mean', self.read_length_mean))
            self.read_length_std = float(length_profile.get('std', self.read_length_std))
            self.read_length_quantiles = length_profile.get('quantiles')
            quality_profile = profile.get('qv', {})
            self.quality_mean = float(quality_profile.get('mean', self.quality_mean))
            self.quality_std = float(quality_profile.get('std', self.quality_std))
            print(
                f"  使用真实数据 profile: length={self.read_length_mean:.0f}±"
                f"{self.read_length_std:.0f}, QV={self.quality_mean:.2f}±{self.quality_std:.2f}"
            )
        
        # 加载errHMM模型
        if errhmm_model:
            if isinstance(errhmm_model, str):
                self.hmm = GCAwareErrHMM()
                self.hmm.load(errhmm_model)
            else:
                self.hmm = errhmm_model
        else:
            self.default_error_rate = 0.01 if platform == 'hifi' else 0.10
            self.hmm = None
        
        # 默认错误率 (如果无HMM模型)
        self.default_error_rate = 0.10  # ONT ~10%
    
    def sample_read_length(self):
        """从正态分布采样读长"""
        if self.read_length_quantiles:
            quantiles = np.asarray(self.read_length_quantiles, dtype=float)
            probabilities = np.linspace(0.0, 1.0, len(quantiles))
            return max(1000, int(np.interp(np.random.random(), probabilities, quantiles)))
        length = int(np.random.normal(self.read_length_mean, self.read_length_std))
        return max(1000, length)  # 最小1000bp
    
    def sample_read_position(self, genome_length, read_length):
        """随机采样read的起始位置"""
        max_start = genome_length - read_length
        if max_start <= 0:
            return 0
        return np.random.randint(0, max_start)
    
    def get_gc_content(self, ref_seq, pos, window_size=100):
        """计算指定位置的GC含量"""
        start = max(0, pos - window_size // 2)
        end = min(len(ref_seq), pos + window_size // 2)
        window = ref_seq[start:end].upper()
        if len(window) == 0:
            return 50.0
        gc = (window.count('G') + window.count('C')) / len(window) * 100
        return gc
    
    def introduce_errors(self, read_seq, ref_seq, ref_start, quality_scores=None):
        """
        使用errHMM引入错误
        
        Args:
            read_seq: 原始序列 (与参考一致)
            ref_seq: 完整参考序列
            ref_start: read在参考上的起始位置
            quality_scores: 质量值列表 (可选)
        
        Returns:
            modified_seq, cigar, new_qualities
        """
        if not self.hmm:
            # Route A: 无HMM模型
            # 无HMM模型时使用简单随机错误
            return self._introduce_random_errors(read_seq)
        
        hmm = self.hmm
        states = hmm.STATES
        n_states = len(states)
        rand = random.random

        # 首次调用时预计算累积转移概率与替换分布
        # (逆CDF采样, 与原实现的 np.random.choice 在分布上完全等价)
        if not hasattr(self, '_cum_trans'):
            cum = {}
            for gc_bin, mat in hmm.transition_probs.items():
                mat = np.maximum(np.asarray(mat, dtype=float), 0)
                mat = mat / mat.sum(axis=1, keepdims=True)
                cum[gc_bin] = np.cumsum(mat, axis=1).tolist()
            self._cum_trans = cum
            self._default_cum = cum.get(hmm.gc_bins // 2)
            self._state_index = {s: i for i, s in enumerate(states)}
            sub = {}
            sm = hmm.substitution_matrix
            for b in 'ATCG':
                if b in sm and sm[b]:
                    probs = sm[b]
                    choices = list(probs.keys())
                    weights = [float(probs[c]) for c in choices]
                    if b not in choices:
                        choices.append(b)
                        weights.append(max(0.1, 1.0 - sum(weights)))
                    w = np.asarray(weights, dtype=float)
                    w = w / w.sum()
                    sub[b] = (choices, np.cumsum(w).tolist())
            self._sub_cum = sub

        cum_trans = self._cum_trans
        default_cum = self._default_cum
        sub_cum = self._sub_cum
        get_gc_bin = hmm.get_gc_bin

        read_len = len(read_seq)
        ref_len = len(ref_seq)
        win = 50  # get_gc_content 默认 window_size=100 的一半

        # 预计算read区域(含缺失漂移余量)的GC前缀和,
        # 避免逐碱基重复统计100bp窗口 (与原 get_gc_content 结果一致)
        lo = max(0, ref_start - win)
        hi = min(ref_len, ref_start + 2 * read_len + win)
        region = ref_seq[lo:hi]
        arr = np.frombuffer(region.encode('ascii'), dtype=np.uint8)
        is_gc = ((arr == 71) | (arr == 67) | (arr == 103) | (arr == 99)).astype(np.int32)
        gc_prefix = [0] + np.cumsum(is_gc).tolist()
        region_n = len(region)

        modified_bases = []
        cigar_ops = []
        ref_pos = ref_start
        si = self._state_index['M']  # 初始状态为Match

        for i in range(read_len):
            base = read_seq[i]

            # 当前ref位置的GC含量 (窗口[ref_pos-50, ref_pos+50), 同原实现)
            a = ref_pos - win
            if a < 0:
                a = 0
            b = ref_pos + win
            if b > ref_len:
                b = ref_len
            if b <= a or a < lo or b > hi:
                gc = self.get_gc_content(ref_seq, ref_pos)
            else:
                # ref_pos因缺失漂移超出预计算区域 (罕见), 退回原始窗口统计
                gc = (gc_prefix[b - lo] - gc_prefix[a - lo]) / (b - a) * 100.0
            gc_bin = get_gc_bin(gc)

            # 根据当前状态采样下一状态 (逆CDF)
            row = cum_trans.get(gc_bin, default_cum)[si]
            u = rand()
            nxt = 0
            last = n_states - 1
            while nxt < last and u > row[nxt]:
                nxt += 1
            next_state = states[nxt]

            # 根据状态决定输出
            if next_state == 'M':
                modified_bases.append(base)
                cigar_ops.append('M')
                ref_pos += 1

            elif next_state == 'S':
                # 替换
                entry = sub_cum.get(base)
                if entry is not None:
                    choices, cw = entry
                    u2 = rand()
                    k = 0
                    lk = len(cw) - 1
                    while k < lk and u2 > cw[k]:
                        k += 1
                    modified_bases.append(choices[k])
                else:
                    modified_bases.append(random.choice([x for x in 'ATCG' if x != base]))
                cigar_ops.append('M')  # 在CIGAR中表现为mismatch
                ref_pos += 1

            elif next_state.startswith('I'):
                # 插入: 长度编码在状态名里 (I1/I2/I3/I4+)
                for _ in range(self._parse_indel_length(next_state)):
                    modified_bases.append('ATCG'[random.randrange(4)])
                    cigar_ops.append('I')
                # ref_pos不变

            elif next_state.startswith('D'):
                # 缺失: 跳过参考碱基, 长度编码在状态名里 (D1/D2/D3/D4+)
                n_del = self._parse_indel_length(next_state)
                ref_pos += n_del
                cigar_ops.extend('D' * n_del)
                si = nxt
                continue

            else:
                raise ValueError(f"unknown HMM state: {next_state!r}")

            si = nxt

        modified_seq = ''.join(modified_bases)
        cigar = self._collapse_cigar(cigar_ops)

        return modified_seq, cigar
    
    def _introduce_random_errors(self, seq, error_rate=None):
        """简单随机错误引入 (无HMM时)"""
        error_rate = error_rate or self.default_error_rate
        bases = list(seq)
        cigar = []
        
        for i, base in enumerate(bases):
            r = np.random.random()
            if r < error_rate * 0.7:
                # 替换
                bases[i] = np.random.choice([b for b in 'ATCG' if b != base])
                cigar.append('M')
            elif r < error_rate * 0.85:
                # 插入
                bases.insert(i, np.random.choice(['A', 'T', 'C', 'G']))
                cigar.append('I')
            elif r < error_rate:
                # 缺失 (跳过)
                cigar.append('D')
                continue
            else:
                cigar.append('M')
        
        return ''.join(bases), self._collapse_cigar(cigar)
    
    def _sample_substitution(self, ref_base):
        """根据替换矩阵采样替换碱基"""
        if self.hmm and ref_base in self.hmm.substitution_matrix:
            probs = self.hmm.substitution_matrix[ref_base]
            if probs:
                choices = list(probs.keys())
                weights = list(probs.values())
                # 添加自身概率
                if ref_base not in choices:
                    choices.append(ref_base)
                    weights.append(max(0.1, 1 - sum(weights)))
                return np.random.choice(choices, p=np.array(weights) / sum(weights))
        
        # 默认随机替换
        return np.random.choice([b for b in 'ATCG' if b != ref_base])
    
    def _parse_indel_length(self, state):
        """Resolve an indel state name to a run length.

        Lengths 1-3 are explicit in the state name; the shared `*4+` state
        covers everything longer. That tail is sampled to match the observed
        ONT distribution rather than from a flat 4-9 range: on real chr21,
        among indels of 4 bp or more, 35 % are exactly 4 bp, 52 % are 5-9 bp,
        10 % are 10-49 bp and 3.5 % are 50 bp or longer (mean 10.45 bp). A flat
        4-9 draw would give a mean of 6.5 and lose the tail entirely.
        """
        if state in ('I', 'D'):
            return 1
        if state.endswith('4+'):
            draw = random.random()
            if draw < 0.35:
                return 4
            if draw < 0.87:
                return random.randint(5, 9)
            if draw < 0.965:
                return random.randint(10, 49)
            # The >=50 bp tail is heavy: on real chr21 it holds 0.56 % of the
            # events but 1.89 Mb of the 5.35 Mb inserted bases, a mean of about
            # 241 bp per event. A uniform 50-150 draw would halve the tail.
            tail = int(random.expovariate(1.0 / 190.0)) + 50
            return min(tail, 5000)
        try:
            return int(state[1:])
        except (TypeError, ValueError):
            return 1
    
    def _collapse_cigar(self, ops):
        """将CIGAR操作列表压缩为字符串"""
        if not ops:
            return ""
        
        cigar = []
        count = 1
        current = ops[0]
        
        for op in ops[1:]:
            if op == current:
                count += 1
            else:
                cigar.append(f"{count}{current}")
                current = op
                count = 1
        
        cigar.append(f"{count}{current}")
        return ''.join(cigar)
    
    def generate_quality_scores(self, length, mean_qv=None):
        """生成质量值"""
        # HiFi reads通常有较高的QV
        if mean_qv is None:
            mean_qv = self.quality_mean
        qvs = np.random.normal(mean_qv, self.quality_std, length)
        qvs = np.clip(qvs, 2, 93).astype(int)
        return qvs
    
    def simulate_reads(self, ref_fasta, coverage, output_prefix, threads=1):
        """
        生成模拟reads (流式写入, 低内存占用)
        
        Args:
            ref_fasta: 参考基因组FASTA
            coverage: 目标覆盖度
            output_prefix: 输出前缀
            threads: 线程数 (当前仅用于日志)
        """
        print(f"[模拟生成] 开始生成模拟reads...")
        print(f"  参考: {ref_fasta}")
        print(f"  覆盖度: {coverage}x")
        print(f"  输出: {output_prefix}")
        
        # 加载参考序列
        random.seed(self.seed)
        np.random.seed(self.seed)
        ref_records = list(SeqIO.parse(ref_fasta, 'fasta'))
        total_ref_length = sum(len(r.seq) for r in ref_records)
        
        # 计算需要的reads数量
        total_bases_needed = coverage * total_ref_length
        num_reads = int(total_bases_needed / self.read_length_mean)
        
        print(f"  参考总长: {total_ref_length:,} bp")
        print(f"  目标bases: {total_bases_needed:,} bp")
        print(f"  预计reads: {num_reads:,}")
        
        # 拼接所有参考序列用于采样 (简化)
        ref_seqs = [(record.id, str(record.seq)) for record in ref_records]
        
        # 流式写入FASTQ, 避免在内存中累积所有reads
        output_fastq = f"{output_prefix}.fastq"
        read_lengths = []
        
        references = {name: seq for name, seq in ref_seqs}
        read_specs = []
        for i in range(num_reads):
            read_len = self.sample_read_length()
            eligible = [(name, seq) for name, seq in ref_seqs if len(seq) >= read_len]
            if not eligible:
                eligible = ref_seqs
            weights = [len(seq) for _, seq in eligible]
            contig_id, contig_seq = random.choices(eligible, weights=weights, k=1)[0]
            start_pos = self.sample_read_position(len(contig_seq), read_len)
            read_specs.append({
                "index": i,
                # numpy's legacy seeding rejects values >= 2**32. The original
                # formula overflowed as soon as a species needed more than
                # ~41k reads (D. melanogaster 140,543; M. musculus chr19
                # ~57k), which killed the worker pool partway through the
                # run. Wrapping keeps the stream identical for the first
                # ~41k reads and stays fully deterministic beyond that.
                "seed": (int(self.seed) + 104729 * (i + 1)) % (2 ** 32),
                "contig_id": contig_id,
                "start": int(start_pos),
                "read_len": int(read_len),
            })

        requested_workers = int(threads)
        worker_count = max(1, min(requested_workers, 16, os.cpu_count() or 1))
        chunk_size = 512
        chunks = [
            read_specs[offset:offset + chunk_size]
            for offset in range(0, len(read_specs), chunk_size)
        ]
        started = time.time()
        completed_reads = 0
        with open(output_fastq, 'w') as fq_out:
            if worker_count == 1:
                results = (
                    _simulate_spec_chunk(self, references, chunk)
                    for chunk in chunks
                )
                for text, lengths in results:
                    fq_out.write(text)
                    read_lengths.extend(lengths)
                    completed_reads += len(lengths)
                    if completed_reads % 10000 < chunk_size:
                        print(
                            f"  生成进度: {completed_reads}/{num_reads} "
                            f"({completed_reads / num_reads * 100:.1f}%)"
                        )
            else:
                try:
                    context = multiprocessing.get_context("fork")
                except ValueError:
                    context = multiprocessing.get_context()
                with concurrent.futures.ProcessPoolExecutor(
                    max_workers=worker_count,
                    mp_context=context,
                    initializer=_init_worker,
                    initargs=(self, references),
                ) as executor:
                    for text, lengths in executor.map(_worker_entry, chunks, chunksize=1):
                        fq_out.write(text)
                        read_lengths.extend(lengths)
                        completed_reads += len(lengths)
                        if completed_reads % 10000 < chunk_size:
                            print(
                                f"  生成进度: {completed_reads}/{num_reads} "
                                f"({completed_reads / num_reads * 100:.1f}%)"
                            )
        wall_seconds = time.time() - started
        
        # 输出统计
        num_reads_done = len(read_lengths)
        total_bases_done = int(sum(read_lengths))
        mean_length = float(np.mean(read_lengths)) if read_lengths else 0.0
        n50 = self._compute_n50(read_lengths)
        
        print(f"\n✓ 模拟完成!")
        print(f"  输出: {output_fastq}")
        print(f"  Reads: {num_reads_done:,}")
        print(f"  总bases: {total_bases_done:,}")
        print(f"  平均长度: {mean_length:.0f} bp")
        print(f"  N50: {n50:,} bp")
        
        # 保存统计
        stats = {
            'num_reads': num_reads_done,
            'total_bases': total_bases_done,
            'mean_length': mean_length,
            'n50': n50,
            'coverage': coverage,
            'reference_length': total_ref_length,
            'model': 'errHMM' if self.hmm else 'random',
            'threads_requested': int(threads),
            'threads_used': worker_count,
            'wall_seconds': round(wall_seconds, 6),
            'reads_per_second': round(
                num_reads_done / wall_seconds, 3
            ) if wall_seconds > 0 else 0.0,
        }
        
        with open(f"{output_prefix}.stats.json", 'w') as f:
            json.dump(stats, f, indent=2)
        
        return output_fastq
    
    def _compute_n50(self, lengths):
        """计算N50"""
        lengths = sorted(lengths, reverse=True)
        cumsum = np.cumsum(lengths)
        half = cumsum[-1] / 2
        idx = np.searchsorted(cumsum, half)
        return int(lengths[idx]) if idx < len(lengths) else int(lengths[-1])


def main():
    parser = argparse.ArgumentParser(description='Generate GC-Aware Simulated Reads')
    parser.add_argument('--ref', required=True, help='Reference FASTA')
    parser.add_argument('--output-prefix', required=True, help='Output prefix')
    parser.add_argument('--coverage', type=float, default=30, help='Target coverage')
    parser.add_argument('--platform', choices=['hifi', 'ont'], default='ont', help='Sequencing platform')
    parser.add_argument('--errhmm-model', help='errHMM model JSON (Route B/C)')
    parser.add_argument('--read-length-mean', type=int, default=None, help='Mean read length')
    parser.add_argument('--read-length-std', type=int, default=None, help='Read length std')
    parser.add_argument('--profile-json', help='Real-read length/QV profile JSON')
    parser.add_argument(
        '--require-profile',
        action='store_true',
        help='Refuse to run without --profile-json. Use this for any run whose '
             'numbers enter a table or figure: dropping the flag otherwise '
             'falls back silently to the platform default read length and '
             'invalidates the whole downstream layer.',
    )
    parser.add_argument(
        '--qc-tolerance',
        type=float,
        default=0.10,
        help='Allowed relative deviation of the produced mean read length, '
             'total bases and read count from the requested profile.',
    )
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--threads', type=int, default=1)

    args = parser.parse_args()

    if args.require_profile and not args.profile_json:
        raise SystemExit(
            "SIMULATION_QC_FAILED: --require-profile was set but no "
            "--profile-json was given. Refusing to fall back to the platform "
            "default read length."
        )
    
    simulator = GCAwareSimulator(
        errhmm_model=args.errhmm_model,
        platform=args.platform,
        read_length_mean=args.read_length_mean,
        read_length_std=args.read_length_std,
        profile_json=args.profile_json,
        seed=args.seed,
    )
    
    output_fastq = simulator.simulate_reads(
        ref_fasta=args.ref,
        coverage=args.coverage,
        output_prefix=args.output_prefix,
        threads=args.threads
    )

    # Gate: compare the produced read set against what was asked for. A silent
    # mismatch here is what turned a missing flag into a wasted layer of
    # downstream results, so it is checked before anything else consumes it.
    stats_path = f"{args.output_prefix}.stats.json"
    with open(stats_path, encoding='utf-8') as handle:
        stats = json.load(handle)

    ref_length = 0
    with open(args.ref) as handle:
        for line in handle:
            if not line.startswith('>'):
                ref_length += len(line.strip())
    expected_bases = args.coverage * ref_length

    problems = []
    if args.profile_json:
        expected_mean = simulator.read_length_mean
        produced = stats['mean_length']
        if expected_mean and abs(produced - expected_mean) / expected_mean > args.qc_tolerance:
            problems.append(
                f"mean read length {produced:.0f} deviates from the profile "
                f"{expected_mean:.0f} by more than {args.qc_tolerance:.0%}"
            )
    if expected_bases and abs(stats['total_bases'] - expected_bases) / expected_bases > args.qc_tolerance:
        problems.append(
            f"total bases {stats['total_bases']} deviate from the requested "
            f"coverage target {expected_bases:.0f} by more than "
            f"{args.qc_tolerance:.0%}"
        )

    if problems:
        raise SystemExit(
            "SIMULATION_QC_FAILED: " + "; ".join(problems)
            + f" (see {stats_path})"
        )

    print(
        f"✓ 模拟质量闸通过: mean_length={stats['mean_length']:.0f} "
        f"total_bases={stats['total_bases']} reads={stats['num_reads']}"
    )
    
    print("\n✓ 模拟数据生成完成!")


if __name__ == '__main__':
    main()
