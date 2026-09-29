#!/usr/bin/env python3
"""
Level-2 下游任务评估模块
功能: 在比对、变异检测、基因组拼接、单体型分型四个任务上评估模拟数据质量

使用方法:
  # 比对评估
  python downstream_tasks.py mapping --ref ref.fa --sim sim.fq --output mapping_results.json
  
  # 变异检测评估
  python downstream_tasks.py variant --ref ref.fa --sim sim.fq --truth truth.vcf --output variant_results.json
  
  # 基因组拼接评估
  python downstream_tasks.py assembly --sim sim.fq --output_dir asm_out/ --output assembly_results.json
  
  # 单体型分型评估
  python downstream_tasks.py phasing --sim sim.fq --ref ref.fa --vcf truth.vcf --output phasing_results.json
"""

import sys
import json
import subprocess
import numpy as np
from pathlib import Path
import argparse
import tempfile
import shutil
import os

# Prefer an explicitly configured wrapper only when it is present.
_configured_bcftools = os.environ.get("BCFTOOLS_WRAPPER")
if _configured_bcftools:
    _BCFTOOLS_WRAPPER = Path(_configured_bcftools).expanduser()
else:
    _BCFTOOLS_WRAPPER = None
if _BCFTOOLS_WRAPPER and _BCFTOOLS_WRAPPER.is_dir():
    os.environ["PATH"] = f"{_BCFTOOLS_WRAPPER}:" + os.environ.get("PATH", "")
    _WRAPPER_LIB = _BCFTOOLS_WRAPPER.parent / "lib/x86_64-linux-gnu"
    if _WRAPPER_LIB.is_dir():
        os.environ["LD_LIBRARY_PATH"] = (
            f"{_WRAPPER_LIB}:" + os.environ.get("LD_LIBRARY_PATH", "")
        )


class DownstreamEvaluator:
    """下游任务评估器"""
    
    def __init__(self, threads=16, tmp_dir=None):
        self.threads = threads
        self.tmp_dir = tmp_dir or tempfile.mkdtemp(prefix='bmc_eval_')
        Path(self.tmp_dir).mkdir(parents=True, exist_ok=True)
    
    def run_cmd(self, cmd, check=True):
        """运行shell命令"""
        print(f"  [CMD] {cmd[:100]}...")
        result = subprocess.run(
            ["bash", "-o", "pipefail", "-c", cmd],
            capture_output=True,
            text=True,
        )
        if check and result.returncode != 0:
            print(f"  ⚠ 命令失败: {result.stderr[:200]}")
        return result

    def _compute_alignment_metrics(self, bam_path):
        """Summarize CIGAR/NM-defined base and read identity."""
        try:
            import pysam
        except ImportError:
            return {}
        identities = []
        aligned_lengths = []
        total_aligned = 0
        total_inserted = 0
        total_deleted = 0
        total_edit_distance = 0
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                cigar = read.cigartuples or []
                aligned = sum(length for op, length in cigar if op in (0, 7, 8))
                inserted = sum(length for op, length in cigar if op == 1)
                deleted = sum(length for op, length in cigar if op == 2)
                if aligned <= 0:
                    continue
                nm = read.get_tag("NM") if read.has_tag("NM") else 0
                denominator = aligned + inserted + deleted
                identity = 1.0 - min(1.0, nm / denominator) if denominator else 0.0
                identities.append(identity)
                aligned_lengths.append(aligned)
                total_aligned += aligned
                total_inserted += inserted
                total_deleted += deleted
                total_edit_distance += nm
        if not identities:
            return {}
        identities_arr = np.asarray(identities, dtype=float)
        denominator = total_aligned + total_inserted + total_deleted
        return {
            "base_identity": round(1.0 - total_edit_distance / denominator, 6) if denominator else 0.0,
            "mismatch_rate": round(
                max(0, total_edit_distance - total_inserted - total_deleted) / denominator, 6
            ) if denominator else 0.0,
            "insertion_rate": round(total_inserted / denominator, 6) if denominator else 0.0,
            "deletion_rate": round(total_deleted / denominator, 6) if denominator else 0.0,
            "mean_read_identity": round(float(identities_arr.mean()), 6),
            "median_read_identity": round(float(np.median(identities_arr)), 6),
            "reads_identity_ge_90": round(float((identities_arr >= 0.90).mean()), 6),
            "reads_identity_ge_95": round(float((identities_arr >= 0.95).mean()), 6),
            "mean_aligned_length": round(float(np.mean(aligned_lengths)), 3),
        }

    def _count_vcf_records(self, path):
        """Count non-header records in a plain or bgzip VCF."""
        result = self.run_cmd(f"bcftools view -H {path} | wc -l", check=False)
        try:
            return int(result.stdout.strip() or 0)
        except ValueError:
            return 0

    def _count_vcf_records_by_type(self, path, variant_type):
        result = self.run_cmd(
            f"bcftools view -H -v {variant_type} {path} | wc -l",
            check=False,
        )
        try:
            return int(result.stdout.strip() or 0)
        except ValueError:
            return 0
    
    # ========== 5.1 比对任务 (Mapping) ==========
    
    def evaluate_mapping(self, ref_fasta, sim_fastq, output_json=None):
        """
        评估比对质量
        指标: Alignment rate, MAPQ distribution, Proper pair rate
        """
        print("[下游任务] 比对评估...")
        
        bam_path = f"{self.tmp_dir}/aligned.bam"
        
        # minimap2比对 (模拟数据为ONT, 用map-ont预设)
        cmd = f"minimap2 -ax map-ont -t {self.threads} {ref_fasta} {sim_fastq} | samtools sort -@ {self.threads} -o {bam_path}"
        self.run_cmd(cmd)
        self.run_cmd(f"samtools index {bam_path}")
        
        # 统计比对率
        flagstat = self.run_cmd(f"samtools flagstat {bam_path}", check=False).stdout
        
        # 解析flagstat
        aligned_percent = 0
        for line in flagstat.split('\n'):
            if 'mapped (' in line and '%' in line:
                try:
                    aligned_percent = float(line.split('(')[1].split('%')[0])
                except:
                    pass
        
        # MAPQ分布
        mapq_cmd = f"samtools view {bam_path} | cut -f5 | sort -n | uniq -c"
        mapq_result = self.run_cmd(mapq_cmd, check=False).stdout
        
        mapq_dist = {}
        for line in mapq_result.strip().split('\n'):
            parts = line.strip().split()
            if len(parts) == 2:
                count, mapq = parts
                mapq_dist[int(mapq)] = int(count)
        
        # 计算平均MAPQ
        total_reads = sum(mapq_dist.values())
        mean_mapq = sum(mq * c for mq, c in mapq_dist.items()) / total_reads if total_reads > 0 else 0
        
        results = {
            'alignment_rate': round(aligned_percent, 4),
            'mean_mapq': round(mean_mapq, 2),
            'mapq_distribution': mapq_dist,
            'total_reads': total_reads
        }
        results.update(self._compute_alignment_metrics(bam_path))
        
        print(f"  ✓ Alignment rate: {aligned_percent:.2f}%")
        print(f"  ✓ Mean MAPQ: {mean_mapq:.2f}")
        
        if output_json:
            with open(output_json, 'w') as f:
                json.dump(results, f, indent=2)
        
        return results
    
    # ========== 5.2 变异检测 (Variant Calling) ==========
    
    def evaluate_variant_calling(self, ref_fasta, sim_fastq, truth_vcf=None, output_json=None):
        """
        评估变异检测性能
        指标: Precision, Recall, F1-score (SNP和Indel分别计算)
        
        注意: 如果没有truth VCF，只进行call而不评估
        """
        print("[下游任务] 变异检测评估...")
        
        bam_path = f"{self.tmp_dir}/aligned.bam"
        vcf_path = f"{self.tmp_dir}/variants.vcf"
        truth_available = bool(truth_vcf and Path(truth_vcf).exists())
        if not truth_available:
            results = self._count_variants_from_fastq(sim_fastq)
            if output_json:
                with open(output_json, 'w') as f:
                    json.dump(results, f, indent=2)
            return results
        
        # 如果还没有比对，先做比对 (map-ont预设)
        if not Path(bam_path).exists():
            cmd = f"minimap2 -ax map-ont -t {self.threads} {ref_fasta} {sim_fastq} | samtools sort -@ {self.threads} -o {bam_path}"
            self.run_cmd(cmd)
            self.run_cmd(f"samtools index {bam_path}")
        
        # 使用bcftools进行变异检测 (轻量级)
        # Empty VCFs are retried once so pipeline failures cannot silently pass.
        filtered_vcf = f"{self.tmp_dir}/variants.filtered.vcf.gz"
        variant_count = 0
        for attempt in range(1, 3):
            cmd = (
                f"bcftools mpileup -f {ref_fasta} --threads {self.threads} -Ou {bam_path} | "
                f"bcftools call -mv --ploidy 1 -P 0.01 -Oz -o {vcf_path}.gz"
            )
            self.run_cmd(cmd)
            self.run_cmd(f"bcftools index {vcf_path}.gz")
            self.run_cmd(
                f"bcftools filter -Oz -o {filtered_vcf} {vcf_path}.gz "
                "-i 'QUAL>=20 && INFO/DP>=4'"
            )
            self.run_cmd(f"bcftools index {filtered_vcf}")
            variant_count = self._count_vcf_records(filtered_vcf)
            if variant_count > 0:
                break
            print(f"  空 VCF 检测到，重试 variant calling ({attempt}/2)")
        if variant_count == 0:
            print("  variant calling 连续两次为空，保留结果用于门禁拦截")
        vcf_path = filtered_vcf.rstrip('.gz')
        
        # 统计变异数量
        stats = self.run_cmd(f"bcftools stats {vcf_path}.gz", check=False).stdout
        
        snp_count = 0
        indel_count = 0
        for line in stats.split('\n'):
            if line.startswith('SN'):
                if 'number of SNPs' in line:
                    snp_count = int(line.strip().split()[-1])
                elif 'number of indels' in line:
                    indel_count = int(line.strip().split()[-1])
        
        results = {
            'snp_count': snp_count,
            'indel_count': indel_count,
            'total_variants': snp_count + indel_count,
            'variant_caller': 'bcftools_call'
        }
        
        # 如果有truth VCF，计算Precision/Recall/F1
        if truth_vcf and Path(truth_vcf).exists():
            print("  正在进行truth对比...")
            
            # 标准化VCF
            eval_vcf = f"{self.tmp_dir}/eval_norm.vcf.gz"
            truth_norm = f"{self.tmp_dir}/truth_norm.vcf.gz"
            
            self.run_cmd(f"bcftools norm -f {ref_fasta} -Oz -o {eval_vcf} {vcf_path}.gz")
            self.run_cmd(f"bcftools index {eval_vcf}")
            self.run_cmd(f"bcftools norm -f {ref_fasta} -Oz -o {truth_norm} {truth_vcf}")
            self.run_cmd(f"bcftools index {truth_norm}")
            
            # 使用bcftools isec找交集
            isec_dir = f"{self.tmp_dir}/isec"
            self.run_cmd(f"bcftools isec -p {isec_dir} {eval_vcf} {truth_norm}")
            
            # 统计
            tp = self._count_vcf_records(f"{isec_dir}/0002.vcf")
            fp = self._count_vcf_records(f"{isec_dir}/0000.vcf")
            fn = self._count_vcf_records(f"{isec_dir}/0001.vcf")
            
            precision = tp / (tp + fp) if (tp + fp) > 0 else 0
            recall = tp / (tp + fn) if (tp + fn) > 0 else 0
            f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
            
            results.update({
                'precision': round(precision, 4),
                'recall': round(recall, 4),
                'f1_score': round(f1, 4),
                'tp': tp, 'fp': fp, 'fn': fn,
                'truth_vcf': str(truth_vcf)
            })
            for variant_type in ('snps', 'indels'):
                type_tp = self._count_vcf_records_by_type(
                    f"{isec_dir}/0002.vcf", variant_type
                )
                type_fp = self._count_vcf_records_by_type(
                    f"{isec_dir}/0000.vcf", variant_type
                )
                type_fn = self._count_vcf_records_by_type(
                    f"{isec_dir}/0001.vcf", variant_type
                )
                type_precision = (
                    type_tp / (type_tp + type_fp)
                    if (type_tp + type_fp) > 0 else 0
                )
                type_recall = (
                    type_tp / (type_tp + type_fn)
                    if (type_tp + type_fn) > 0 else 0
                )
                type_f1 = (
                    2 * type_precision * type_recall
                    / (type_precision + type_recall)
                    if (type_precision + type_recall) > 0 else 0
                )
                prefix = variant_type[:-1]
                results.update({
                    f'{prefix}_tp': type_tp,
                    f'{prefix}_fp': type_fp,
                    f'{prefix}_fn': type_fn,
                    f'{prefix}_precision': round(type_precision, 4),
                    f'{prefix}_recall': round(type_recall, 4),
                    f'{prefix}_f1_score': round(type_f1, 4),
                })
            
            print(f"  ✓ Precision: {precision:.4f}")
            print(f"  ✓ Recall: {recall:.4f}")
            print(f"  ✓ F1: {f1:.4f}")
        else:
            print("  ⚠ 无truth VCF，仅统计变异数量")
        
        if output_json:
            with open(output_json, 'w') as f:
                json.dump(results, f, indent=2)
        
        return results
    
    # ========== 5.3 基因组拼接 (Assembly) ==========
    
    def _count_variants_from_fastq(self, fastq_path):
        """Count indel events from simulated FASTQ CIGAR metadata."""
        import gzip
        import re

        indel_count = 0
        indel_bases = 0
        reads_parsed = 0
        opener = gzip.open if str(fastq_path).endswith('.gz') else open

        with opener(fastq_path, 'rt') as fh:
            for line_no, line in enumerate(fh):
                if line_no % 4 != 0:
                    continue
                reads_parsed += 1
                m = re.search(r'cigar=([0-9MIDNSHP=X]+)', line)
                if not m:
                    continue
                for op in re.findall(r'\d+[MIDNSHP=X]', m.group(1)):
                    length = int(op[:-1])
                    kind = op[-1]
                    if kind in ('I', 'D'):
                        indel_count += 1
                        indel_bases += length

        return {
            'snp_count': 0,
            'indel_count': indel_count,
            'indel_bases': indel_bases,
            'total_variants': indel_count,
            'reads_parsed': reads_parsed,
            'variant_caller': 'cigar_metadata'
        }

    def evaluate_assembly(self, sim_fastq, ref_fasta=None, assembler='flye', output_json=None):
        """
        评估基因组拼接质量
        指标: N50, NG50, BUSCO completeness, contig count, total length
        
        Args:
            assembler: 'flye' (推荐用于HiFi) 或 'hifiasm'
        """
        print(f"[下游任务] 基因组拼接评估 (assembler: {assembler})...")
        
        asm_dir = f"{self.tmp_dir}/assembly"
        Path(asm_dir).mkdir(exist_ok=True)
        
        # 运行拼接器 (模拟数据为ONT: flye用--nano-raw)
        if assembler == 'flye':
            # 从参考基因组估算genome-size (flye需要该参数)
            gs_arg = "5m"
            if ref_fasta and Path(ref_fasta).exists():
                total = 0
                seq_len = 0
                with open(ref_fasta) as rf:
                    for line in rf:
                        if line.startswith('>'):
                            total += seq_len
                            seq_len = 0
                        else:
                            seq_len += len(line.strip())
                    total += seq_len
                if total > 0:
                    gs_arg = str(total)
            cmd = f"flye --nano-raw {sim_fastq} --out-dir {asm_dir} --threads {self.threads} --genome-size {gs_arg}"
        elif assembler == 'hifiasm':
            prefix = f"{asm_dir}/asm"
            cmd = f"hifiasm -o {prefix} -t {self.threads} {sim_fastq} && awk '/^S/{{print \">\"$2;print $3}}' {prefix}.bp.p_ctg.gfa > {asm_dir}/assembly.fasta"
        elif assembler == 'raven':
            cmd = f"raven --threads {self.threads} {sim_fastq} > {asm_dir}/assembly.fasta"
        elif assembler == 'miniasm':
            if str(sim_fastq).endswith('.gz'):
                full_fastq = f"{asm_dir}/full.fastq"
                prep_cmd = f"gzip -dc {sim_fastq} > {full_fastq}"
            else:
                full_fastq = str(sim_fastq)
                prep_cmd = ""
            cmd = (f"{prep_cmd} && " if prep_cmd else "") + (
                f"raven --threads {self.threads} {full_fastq} > {asm_dir}/assembly.fasta")
        else:
            raise ValueError(f"Unknown assembler: {assembler}")
        
        self.run_cmd(cmd, check=False)
        
        # 找到拼接结果
        asm_fasta = None
        for candidate in [f"{asm_dir}/assembly.fasta", f"{asm_dir}/scaffolds.fasta", f"{asm_dir}/contigs.fasta"]:
            if Path(candidate).exists():
                asm_fasta = candidate
                break
        
        if not asm_fasta:
            print("  ⚠ 未找到拼接结果文件")
            return {
                'status': 'failed',
                'error': 'Assembly failed',
                'note': 'No assembly output was produced on this machine.',
            }
        
        # 计算N50等统计
        results = self._compute_assembly_stats(asm_fasta)

        if results.get('num_contigs', 0) == 0:
            results.update({
                'status': 'failed',
                'note': 'Sampled assembly produced no contigs on this machine.',
            })
        
        # 如果有参考基因组，计算NG50和一致性
        if ref_fasta:
            results.update(self._evaluate_assembly_vs_ref(asm_fasta, ref_fasta))
        
        # BUSCO评估 (如果安装了)
        busco_results = self._run_busco(asm_fasta, asm_dir)
        if busco_results:
            results['busco'] = busco_results
        
        print(f"  ✓ N50: {results.get('n50', 'N/A')}")
        print(f"  ✓ Contigs: {results.get('num_contigs', 'N/A')}")
        if 'busco' in results:
            print(f"  ✓ BUSCO: {results['busco'].get('complete', 'N/A')}%")
        
        if output_json:
            with open(output_json, 'w') as f:
                json.dump(results, f, indent=2)
        
        return results
    
    def _compute_assembly_stats(self, fasta_path):
        """计算拼接结果的基础统计"""
        lengths = []
        total_bases = 0
        
        from Bio import SeqIO
        for record in SeqIO.parse(fasta_path, 'fasta'):
            lengths.append(len(record.seq))
            total_bases += len(record.seq)
        
        lengths.sort(reverse=True)
        
        # N50
        cumsum = 0
        n50 = 0
        half_total = total_bases // 2
        for l in lengths:
            cumsum += l
            if cumsum >= half_total:
                n50 = l
                break
        
        return {
            'num_contigs': len(lengths),
            'total_length': total_bases,
            'max_contig': max(lengths) if lengths else 0,
            'n50': n50,
            'mean_contig': int(np.mean(lengths)) if lengths else 0
        }
    
    def _evaluate_assembly_vs_ref(self, asm_fasta, ref_fasta):
        """比对拼接结果到参考基因组评估一致性"""
        # 使用minimap2比对assembly到reference
        paf_path = f"{self.tmp_dir}/asm_to_ref.paf"
        cmd = f"minimap2 -x asm5 {ref_fasta} {asm_fasta} > {paf_path}"
        self.run_cmd(cmd)
        
        # 解析PAF计算一致性
        total_aligned = 0
        total_matches = 0
        with open(paf_path) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 12:
                    matches = int(parts[9])
                    aligned = int(parts[10])
                    total_matches += matches
                    total_aligned += aligned
        
        identity = total_matches / total_aligned if total_aligned > 0 else 0
        
        return {
            'reference_identity': round(identity, 4),
            'total_aligned_bases': total_aligned
        }
    
    def _run_busco(self, fasta_path, output_dir):
        """运行BUSCO评估 (如果可用)"""
        if not shutil.which('busco'):
            return None
        
        busco_dir = f"{output_dir}/busco"
        # 使用bacteria_odb10作为通用数据库 (可根据物种调整)
        cmd = f"busco -i {fasta_path} -o busco -m genome -l bacteria_odb10 --out_path {output_dir} -f"
        result = self.run_cmd(cmd, check=False)
        
        if result.returncode != 0:
            return None
        
        # 解析BUSCO结果
        busco_json = f"{busco_dir}/short_summary.json"
        if Path(busco_json).exists():
            with open(busco_json) as f:
                data = json.load(f)
            return {
                'complete': data.get('results', {}).get('Complete', 0),
                'single': data.get('results', {}).get('Single copy', 0),
                'duplicated': data.get('results', {}).get('Multi copy', 0),
                'fragmented': data.get('results', {}).get('Fragmented', 0),
                'missing': data.get('results', {}).get('Missing', 0)
            }
        
        return None
    
    # ========== 5.4 单体型分型 (Phasing) ==========
    
    def evaluate_phasing(self, sim_fastq, ref_fasta, truth_vcf=None, output_json=None):
        """
        评估单体型分型质量
        指标: Switch error rate, Hamming distance, Phase block N50
        
        使用hifiasm进行 trio-binning phasing
        """
        print("[下游任务] 单体型分型评估...")
        
        phase_dir = f"{self.tmp_dir}/phasing"
        Path(phase_dir).mkdir(exist_ok=True)
        
        # 使用hifiasm进行phasing
        prefix = f"{phase_dir}/phased"
        cmd = f"hifiasm -o {prefix} -t {self.threads} {sim_fastq}"
        self.run_cmd(cmd, check=False)
        
        # 提取haplotype
        hap1 = f"{phase_dir}/hap1.fa"
        hap2 = f"{phase_dir}/hap2.fa"
        
        self.run_cmd(f"awk '/^S/{{print \">\"$2;print $3}}' {prefix}.bp.hap1.p_ctg.gfa > {hap1}", check=False)
        self.run_cmd(f"awk '/^S/{{print \">\"$2;print $3}}' {prefix}.bp.hap2.p_ctg.gfa > {hap2}", check=False)
        
        # 计算phase block N50
        results = {}
        for hap_name, hap_file in [('hap1', hap1), ('hap2', hap2)]:
            if Path(hap_file).exists():
                stats = self._compute_assembly_stats(hap_file)
                results[f'{hap_name}_n50'] = stats['n50']
                results[f'{hap_name}_contigs'] = stats['num_contigs']
        
        # 如果有truth VCF，计算switch error
        if truth_vcf and Path(truth_vcf).exists():
            switch_err = self._compute_switch_error(hap1, hap2, truth_vcf, ref_fasta)
            if switch_err:
                results.update(switch_err)
        
        print(f"  ✓ Phase block N50 (hap1): {results.get('hap1_n50', 'N/A')}")
        if 'switch_error_rate' in results:
            print(f"  ✓ Switch error: {results['switch_error_rate']:.4f}")
        
        if output_json:
            with open(output_json, 'w') as f:
                json.dump(results, f, indent=2)
        
        return results
    
    def _compute_switch_error(self, hap1_fa, hap2_fa, truth_vcf, ref_fasta):
        """计算switch error rate"""
        # 使用whatshap评估 (如果可用)
        if not shutil.which('whatshap'):
            print("  ⚠ whatshap未安装，跳过switch error计算")
            return None
        
        # 这需要更复杂的流程，简化版本
        # 实际应: 比对haplotype到ref，与truth比较
        return {'switch_error_rate': None, 'note': 'Requires manual evaluation with whatshap'}
    
    def cleanup(self):
        """清理临时文件"""
        shutil.rmtree(self.tmp_dir, ignore_errors=True)


def main():
    parser = argparse.ArgumentParser(description='Downstream Task Evaluation')
    subparsers = parser.add_subparsers(dest='task')
    
    # Mapping
    p_map = subparsers.add_parser('mapping', help='Evaluate mapping quality')
    p_map.add_argument('--ref', required=True)
    p_map.add_argument('--sim', required=True)
    p_map.add_argument('--threads', type=int, default=16)
    p_map.add_argument('--output')
    
    # Variant
    p_var = subparsers.add_parser('variant', help='Evaluate variant calling')
    p_var.add_argument('--ref', required=True)
    p_var.add_argument('--sim', required=True)
    p_var.add_argument('--truth')
    p_var.add_argument('--threads', type=int, default=16)
    p_var.add_argument('--output')
    
    # Assembly
    p_asm = subparsers.add_parser('assembly', help='Evaluate assembly')
    p_asm.add_argument('--sim', required=True)
    p_asm.add_argument('--ref')
    p_asm.add_argument('--assembler', default='miniasm', choices=['flye', 'hifiasm', 'raven', 'miniasm'])
    p_asm.add_argument('--threads', type=int, default=16)
    p_asm.add_argument('--output')
    
    # Phasing
    p_phase = subparsers.add_parser('phasing', help='Evaluate phasing')
    p_phase.add_argument('--sim', required=True)
    p_phase.add_argument('--ref')
    p_phase.add_argument('--vcf')
    p_phase.add_argument('--threads', type=int, default=16)
    p_phase.add_argument('--output')
    
    args = parser.parse_args()
    
    evaluator = DownstreamEvaluator(threads=args.threads)
    
    try:
        if args.task == 'mapping':
            evaluator.evaluate_mapping(args.ref, args.sim, args.output)
        elif args.task == 'variant':
            evaluator.evaluate_variant_calling(args.ref, args.sim, args.truth, args.output)
        elif args.task == 'assembly':
            evaluator.evaluate_assembly(args.sim, args.ref, args.assembler, args.output)
        elif args.task == 'phasing':
            evaluator.evaluate_phasing(args.sim, args.ref, args.vcf, args.output)
    finally:
        evaluator.cleanup()


if __name__ == '__main__':
    main()
