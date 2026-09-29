#!/usr/bin/env python3
"""
errHMM training module
功能: 从真实比对结果中提取GC-aware错误模式，训练HMM转移矩阵

GC-Aware Transition Design:
  - 将参考基因组按100bp窗口分割
  - 计算每个窗口的GC含量，归入GC_bin (0-9, 每10%一个bin)
  - 对每个GC_bin分别统计错误状态转移频率
  - 输出10个转移概率矩阵

状态定义:
  M (Match): 正确匹配
  S (Substitution): 碱基替换 (12种转换)
  I (Insertion): 插入
  D1-D4+: 不同长度的缺失
"""

import sys
import json
import numpy as np
import pysam
from collections import defaultdict, Counter
from Bio import SeqIO
import pickle
import argparse
from pathlib import Path


class GCAwareErrHMM:
    """GC感知的错误HMM模型"""
    
    # 状态定义
    # Indels are encoded with the run length in the state name, for insertions
    # and deletions alike. The original model had a single 'I' state, which made
    # the I row's self-transition unobservable: every insertion was then emitted
    # as exactly 1 bp. Real ONT chr21 has a 3.83 bp mean insertion length
    # (55.6 % 1 bp, 19.8 % 2 bp, ... 0.6 % >=50 bp) and 3.78 % inserted bases,
    # against 1.00 bp and 0.99 % under the old encoding.
    STATES = ['M', 'S', 'I1', 'I2', 'I3', 'I4+', 'D1', 'D2', 'D3', 'D4+']
    
    def __init__(
        self,
        gc_bins=10,
        window_size=100,
        pseudocount=1.0,
        prior_strength=10.0,
    ):
        """
        初始化errHMM
        
        Args:
            gc_bins: GC含量分箱数 (默认10: 0-9%, 10-19%, ..., 90-99%)
            window_size: GC含量计算窗口大小
            pseudocount: 兼容参数；用于计算 pooled prior 的均匀化伪计数
            prior_strength: 每个 GC-bin 向全局状态转移分布收缩的权重
        """
        self.gc_bins = gc_bins
        self.window_size = window_size
        self.pseudocount = pseudocount
        self.prior_strength = prior_strength
        
        # 为每个GC_bin存储转移计数矩阵
        # transition_counts[gc_bin][from_state][to_state] = count
        self.transition_counts = defaultdict(lambda: defaultdict(Counter))
        
        # 训练好的转移概率矩阵
        # transition_probs[gc_bin] = numpy array (n_states x n_states)
        self.transition_probs = {}
        
        # 发射概率 (替换矩阵)
        self.substitution_matrix = defaultdict(Counter)

        # 训练元数据
        self.training_metadata = {}
    
    def get_gc_bin(self, gc_percent):
        """将GC百分比映射到bin索引"""
        bin_idx = int(gc_percent / 100.0 * self.gc_bins)
        return min(bin_idx, self.gc_bins - 1)  # 确保不超过最大值
    
    def compute_window_gc(self, ref_seq, window_start, window_end):
        """计算参考序列窗口的GC含量"""
        window = ref_seq[window_start:window_end].upper()
        if len(window) == 0:
            return 50.0  # 默认中值
        gc_count = window.count('G') + window.count('C')
        return (gc_count / len(window)) * 100.0
    
    def precompute_gc_bins(self, ref_seq):
        """预计算整条参考序列的GC_bin数组 (每条染色体只算一次)"""
        seq_len = len(ref_seq)
        gc_bins_seq = []
        for i in range(0, seq_len, self.window_size):
            gc = self.compute_window_gc(ref_seq, i, min(i + self.window_size, seq_len))
            gc_bins_seq.append(self.get_gc_bin(gc))
        return gc_bins_seq

    def _indel_state(self, kind: str, length: int) -> str:
        """Name an indel state by its run length.

        Runs of 4 or more share one state (I4+/D4+), matching the original
        deletion encoding. If the active state space has been simplified to a
        single generic state for this kind, fall back to it.
        """
        candidate = f"{kind}4+" if length > 3 else f"{kind}{length}"
        if candidate in self.STATES:
            return candidate
        if kind in self.STATES:          # simplified space, e.g. ['M','S','I','D']
            return kind
        raise ValueError(
            f"state {candidate!r} not in model state space {self.STATES!r}"
        )

    def cigar_to_state_sequence(self, read, ref_seq, gc_bins_seq, ref_start):
        """
        将CIGAR字符串转换为状态序列

        Args:
            read: pysam比对记录
            ref_seq: 参考序列
            gc_bins_seq: 预计算的GC_bin数组 (见 precompute_gc_bins)
            ref_start: 比对起始位置

        Returns:
            list of (state, gc_bin, ref_pos) tuples
        """
        states = []
        ref_pos = ref_start
        query_pos = 0
        # query_alignment_sequence is already reverse-complemented for
        # reverse-strand alignments and excludes soft/hard clipping.
        query_seq = read.query_alignment_sequence or ''
        for cigar_op, length in read.cigartuples:
            if cigar_op == 0:  # M (match/mismatch)
                for i in range(length):
                    pos = ref_pos + i
                    gc_bin = gc_bins_seq[pos // self.window_size] if pos // self.window_size < len(gc_bins_seq) else self.gc_bins // 2
                    
                    # 判断是Match还是Substitution
                    ref_base = ref_seq[pos].upper() if pos < len(ref_seq) else 'N'
                    query_base = (
                        query_seq[query_pos + i].upper()
                        if query_pos + i < len(query_seq)
                        else 'N'
                    )
                    
                    if ref_base == query_base:
                        states.append(('M', gc_bin, pos))
                    else:
                        states.append(('S', gc_bin, pos))
                        if ref_base in 'ACGT' and query_base in 'ACGT':
                            self.substitution_matrix[ref_base][query_base] += 1
                
                ref_pos += length
                query_pos += length
                
            elif cigar_op == 1:  # I (insertion)
                gc_bin = gc_bins_seq[ref_pos // self.window_size] if ref_pos // self.window_size < len(gc_bins_seq) else self.gc_bins // 2
                state = self._indel_state('I', length)
                states.append((state, gc_bin, ref_pos))
                query_pos += length

            elif cigar_op == 2:  # D (deletion)
                gc_bin = gc_bins_seq[ref_pos // self.window_size] if ref_pos // self.window_size < len(gc_bins_seq) else self.gc_bins // 2
                state = self._indel_state('D', length)
                states.append((state, gc_bin, ref_pos))
                ref_pos += length
                
            elif cigar_op == 4 or cigar_op == 5:  # S/H (soft/hard clip)
                # query_alignment_sequence excludes clipped bases.
                pass
        
        return states
    
    def train_from_bam(self, bam_path, ref_fasta, min_mapq=20, max_reads=0):
        """
        从BAM文件训练errHMM
        
        Args:
            bam_path: 比对后的BAM文件路径
            ref_fasta: 参考基因组FASTA路径
            min_mapq: 最小比对质量
            max_reads: 最多使用的合格reads数；0 表示不限制
        """
        print(f"[train_errhmm] 开始训练...")
        print(f"  BAM: {bam_path}")
        print(f"  REF: {ref_fasta}")
        print(f"  MAPQ >= {min_mapq}")
        print(f"  GC-bin prior strength: {self.prior_strength}")
        
        # 加载参考序列
        ref_dict = {}
        for record in SeqIO.parse(ref_fasta, 'fasta'):
            ref_dict[record.id] = str(record.seq)
        print(f"  加载了 {len(ref_dict)} 条参考序列")

        # 预计算每条染色体的GC_bin数组 (避免对每条read重复扫描整条染色体)
        gc_bins_dict = {name: self.precompute_gc_bins(seq) for name, seq in ref_dict.items()}
        print(f"  GC_bin数组预计算完成")
        
        # 打开BAM文件
        bam = pysam.AlignmentFile(bam_path, 'rb')
        
        total_reads = 0
        processed_reads = 0
        
        for read in bam:
            total_reads += 1
            
            if read.is_unmapped or read.mapping_quality < min_mapq:
                continue
            if read.is_secondary or read.is_supplementary or read.is_qcfail or read.is_duplicate:
                continue
            if read.query_sequence is None or read.query_length < 100:
                continue
            
            ref_name = read.reference_name
            if ref_name not in ref_dict:
                continue
            
            ref_seq = ref_dict[ref_name]
            states = self.cigar_to_state_sequence(read, ref_seq, gc_bins_dict[ref_name], read.reference_start)
            
            if len(states) < 2:
                continue
            
            # 记录状态转移
            for i in range(len(states) - 1):
                from_state, gc_bin, _ = states[i]
                to_state, _, _ = states[i + 1]
                self.transition_counts[gc_bin][from_state][to_state] += 1
            
            processed_reads += 1

            if max_reads and processed_reads >= max_reads:
                print(f"  达到训练上限: {processed_reads} reads")
                break
            
            if total_reads % 100000 == 0:
                print(f"  处理中... {total_reads} reads, {processed_reads} 合格")
        
        bam.close()
        
        print(f"\n[train_errhmm] 训练完成: {processed_reads}/{total_reads} reads")
        self.training_metadata = {
            'bam': str(bam_path),
            'reference': str(ref_fasta),
            'min_mapq': int(min_mapq),
            'reads_total': int(total_reads),
            'reads_processed': int(processed_reads),
            'max_reads': int(max_reads),
            'gc_bins': int(self.gc_bins),
            'window_size': int(self.window_size),
        }
        
        # 计算转移概率 (带Laplace平滑)
        self._compute_transition_probabilities()
        
        # 计算替换概率
        self._compute_substitution_probs()
    
    def _compute_transition_probabilities(self):
        """使用全局经验分布作为 prior，计算稀疏 GC-bin 的稳健转移矩阵。

        直接对空 GC-bin 做 Laplace(1) 会把整行变成近似均匀分布，模拟器在该
        bin 中会以 6/7 的概率产生错误。这里先用所有 GC-bin 的 pooled counts
        估计全局转移分布，再把每个 bin 向该分布收缩。空 bin 因而退化为全局
        分布，低覆盖 bin 也不会被少量计数主导。
        """
        print("\n[train_errhmm] 计算层次平滑转移概率矩阵...")

        n_states = len(self.STATES)
        global_counts = np.zeros((n_states, n_states), dtype=float)
        raw_counts = {}

        for gc_bin in range(self.gc_bins):
            counts = self.transition_counts[gc_bin]
            count_matrix = np.zeros((n_states, n_states), dtype=float)
            for i, from_state in enumerate(self.STATES):
                for j, to_state in enumerate(self.STATES):
                    count_matrix[i, j] = counts[from_state].get(to_state, 0)
            raw_counts[gc_bin] = count_matrix
            global_counts += count_matrix

        # 每个状态单独归一化；完全没有观测的状态使用 uniform prior。
        global_prior = np.zeros_like(global_counts)
        for i in range(n_states):
            row = global_counts[i]
            if row.sum() > 0:
                global_prior[i] = row / row.sum()
            else:
                global_prior[i] = 1.0 / n_states

        # 正观测尽量保留 GC 条件信号；空 bin 完全回退到全局经验分布。
        for gc_bin in range(self.gc_bins):
            raw = raw_counts[gc_bin]
            # pseudocount 仅用于避免 pooled prior 中某行完全为零；
            # prior_strength 控制池化先验的权重。
            posterior = raw + self.prior_strength * global_prior
            row_sums = posterior.sum(axis=1, keepdims=True)
            self.transition_probs[gc_bin] = posterior / np.where(row_sums == 0, 1, row_sums)

            print(
                f"  GC_bin {gc_bin}: {int(raw.sum())} observed transitions "
                f"(prior weight={self.prior_strength:g})"
            )

        self.training_metadata.update({
            'smoothing': 'hierarchical_empirical_prior',
            'prior_strength': float(self.prior_strength),
            'global_transition_counts': global_counts.tolist(),
        })
        print("✓ 层次平滑转移概率矩阵计算完成")
    
    def _compute_substitution_probs(self):
        """计算碱基替换概率"""
        for ref_base, counts in self.substitution_matrix.items():
            total = sum(counts.values())
            if total > 0:
                for alt_base in counts:
                    self.substitution_matrix[ref_base][alt_base] /= total
    
    def get_transition_matrix(self, gc_percent):
        """获取指定GC含量对应的转移矩阵"""
        gc_bin = self.get_gc_bin(gc_percent)
        return self.transition_probs.get(gc_bin, self.transition_probs.get(self.gc_bins // 2))
    
    def save(self, output_path):
        """保存模型到JSON"""
        model_data = {
            'model_version': 2,
            'gc_bins': self.gc_bins,
            'window_size': self.window_size,
            'pseudocount': self.pseudocount,
            'prior_strength': self.prior_strength,
            'training_metadata': self.training_metadata,
            'states': self.STATES,
            'transition_probs': {
                str(k): v.tolist() for k, v in self.transition_probs.items()
            },
            'substitution_probs': {
                k: dict(v) for k, v in self.substitution_matrix.items()
            }
        }
        
        with open(output_path, 'w') as f:
            json.dump(model_data, f, indent=2)
        
        print(f"✓ 模型已保存: {output_path}")
    
    def load(self, model_path):
        """从JSON加载模型"""
        with open(model_path, 'r') as f:
            data = json.load(f)
        
        self.gc_bins = data['gc_bins']
        self.window_size = data['window_size']
        self.pseudocount = data.get('pseudocount', 1.0)
        self.prior_strength = data.get('prior_strength', 10.0)
        self.training_metadata = data.get('training_metadata', {})
        self.STATES = data['states']
        
        for k, v in data['transition_probs'].items():
            self.transition_probs[int(k)] = np.array(v)
        
        for k, v in data['substitution_probs'].items():
            self.substitution_matrix[k] = Counter(v)
        
        print(f"✓ 模型已加载: {model_path}")
    
    def print_transition_summary(self):
        """打印转移矩阵摘要"""
        print("\n===== errHMM 转移概率摘要 =====")
        for gc_bin in sorted(self.transition_probs.keys()):
            print(f"\nGC_bin {gc_bin} ({gc_bin * 10}-{(gc_bin + 1) * 10}%):")
            matrix = self.transition_probs[gc_bin]
            
            # 打印表头
            header = "      " + " ".join(f"{s:>6}" for s in self.STATES)
            print(header)
            
            for i, from_state in enumerate(self.STATES):
                row_str = f"{from_state:>4}  " + " ".join(f"{matrix[i,j]:>6.4f}" for j in range(len(self.STATES)))
                print(row_str)


def main():
    parser = argparse.ArgumentParser(description='Train GC-Aware errHMM from BAM alignment')
    parser.add_argument('--bam', required=True, help='Input BAM file')
    parser.add_argument('--ref', required=True, help='Reference FASTA')
    parser.add_argument('--output', required=True, help='Output JSON model path')
    parser.add_argument('--gc-bins', type=int, default=10, help='Number of GC bins')
    parser.add_argument('--nongc', action='store_true', help='Collapse all GC bins into one (Route B)')
    parser.add_argument('--window-size', type=int, default=100, help='GC window size')
    parser.add_argument('--min-mapq', type=int, default=20, help='Minimum MAPQ')
    parser.add_argument(
        '--max-reads',
        type=int,
        default=0,
        help='Maximum accepted reads to train on (0 = all)',
    )
    parser.add_argument(
        '--pseudocount',
        type=float,
        default=None,
        help='Deprecated alias for --prior-strength',
    )
    parser.add_argument(
        '--prior-strength',
        type=float,
        default=10.0,
        help='Weight of the pooled empirical transition prior',
    )
    parser.add_argument('--summary', action='store_true', help='Print transition matrix summary')
    
    args = parser.parse_args()
    
    # 训练模型
    gc = 1 if args.nongc else args.gc_bins
    prior_strength = args.prior_strength
    if args.pseudocount is not None:
        prior_strength = args.pseudocount
        print('[WARN] --pseudocount is deprecated; use --prior-strength')

    hmm = GCAwareErrHMM(
        gc_bins=gc,
        window_size=args.window_size,
        pseudocount=1.0,
        prior_strength=prior_strength,
    )
    hmm.train_from_bam(
        args.bam,
        args.ref,
        min_mapq=args.min_mapq,
        max_reads=args.max_reads,
    )
    
    # 保存模型
    hmm.save(args.output)
    
    # 打印摘要
    if args.summary:
        hmm.print_transition_summary()
    
    print("\n✓ 训练完成!")


if __name__ == '__main__':
    main()
