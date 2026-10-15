#!/usr/bin/env python3
"""Write the W3.8 Chinese evaluation response as Markdown."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def format_p(value: float) -> str:
    value = float(value)
    return "<0.0001" if value < 0.0001 else f"{value:.4f}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--stats", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    stats = json.loads(args.stats.read_text(encoding="utf-8"))
    second = stats["second_cohort"]
    pooled = stats["all_seed_sensitivity"]
    supported = bool(second["supported"])
    status = (
        "第二轮独立队列达到预设终点，支持把拟南芥中的 GC 条件化效应写成"
        "单物种确认性结果。"
        if supported
        else
        "第二轮独立队列没有达到预设终点，因此我们仍按未确认结果报告，"
        "没有事后修改阈值。"
    )

    lines = [
        "# W3.8 对评估意见的回应",
        "",
        "这轮修改没有停在措辞层面。我们重排了 Results 的论证顺序，补充了"
        "模型方法、计算性能和仓库说明，并新增了一个在运行前锁定终点与分析"
        "规则的拟南芥第二轮独立复现队列。原来的六物种 gate 和第一轮复现"
        "结果均未追溯修改。",
        "",
        "## 摘要",
        "",
        "Background 现在按“错误依赖局部 GC -> 常用模拟器主要暴露全局误差"
        "-> 常用指标不一定测量局部 GC 误差”展开。Results 压缩了实验细节，"
        "保留受控消融、六物种范围、独立复现和下游一致性。Conclusions 与"
        "Background 提出的两个问题逐项对应。",
        "",
        "## 背景",
        "",
        "补充了 GC 影响长读长错误的机制解释，包括局部碱基组成、孔道信号与"
        "动力学、二级结构、重复区域比对歧义。Figure 1 现在有明确的分析目的、"
        "指标含义和主要观察。对 context-aware simulation 的比较已补齐："
        "明确说明 CycSim 是上下文感知长读长模拟器，并区分其面向 mapping "
        "参数优化的目标与本文的 GC-bin 条件化状态转移和 GC 分层评估。",
        "",
        "## 结果",
        "",
        "Results 已按证据顺序重排：模型与测量问题、受控 1-bin 对 GC-aware"
        "消融、第二轮独立复现、六物种范围和预注册 gate、不确定性、独立 indel"
        "消融、外部工具与下游一致性、计算性能、结论边界。每一块都说明它能"
        "支持什么以及不能支持什么。新增 Table 6 汇总 discovery、第一轮复现、"
        "第二轮复现和 12-seed 敏感性分析。",
        "",
        "## 方法",
        "",
        "补充了 `M/S/I1-I3/I4+/D1-D3/D4+` 状态定义、transition matrix 的"
        "含义、训练计数与 hierarchical smoothing、读段生成流程和 worker"
        "并行机制。软件版本移到独立的 `Software and versions` 小节，不再插入"
        "核心方法叙述中间。",
        "",
        "## 代码仓库",
        "",
        "README 增加了 reviewer quick start、核心模块关系和计算性能命令。"
        "`--workers` 是正式参数，`--threads` 只保留为兼容别名。新增离线小"
        "fixture、参数拼接回归测试、资源 benchmark 和第二轮复现统计脚本。"
        "`--quick-reproduce` 的数据准备参数已修复，并加了测试，避免再次出现"
        "把续行参数当成独立命令执行的运行时错误。",
        "",
        "## 新增的独立复现",
        "",
        "协议在运行前固定：A. thaliana、10x、profile-matched、1-bin 对 5-bin、"
        "新 seeds 77、88、99、111、122、133，主要终点为 seed 内 5-bin 减"
        "1-bin 的 GC-curve Pearson r，双侧 paired t-test，alpha = 0.05，"
        "一次性分析，不根据结果改 seed、终点或 bin 数。",
        "",
        f"结果：mean delta r = {float(second['mean_delta']):.4f}，95% CI "
        f"{float(second['ci95_low']):.4f} 到 {float(second['ci95_high']):.4f}，"
        f"p = {float(second['paired_t_p']):.4f}。{status}",
        "",
        f"把 discovery 3 个 seed、第一轮 3 个 seed 和第二轮 6 个 seed 合并后，"
        f"12-seed 敏感性分析为 mean delta r = {float(pooled['mean_delta']):.4f}，"
        f"95% CI {float(pooled['ci95_low']):.4f} 到 "
        f"{float(pooled['ci95_high']):.4f}，p "
        f"{format_p(pooled['paired_t_p'])}。合并分析只作为敏感性分析，"
        "不替代原预注册 gate。",
        "",
        "## 如何回应 CycSim 意见",
        "",
        "CycSim 是上下文感知长读长模拟器，和我们有直接关系。新稿在 Background、"
        "Discussion 和 Table 10 中明确写出：CycSim 学习序列上下文依赖的错误谱，"
        "主要服务于 mapping 参数优化；GCerrHMM 条件化的是 GC-bin 的 error-state "
        "transition，并用 GC 分层曲线和预注册规则评价。由于本文的跨工具面板已在"
        "统一 coverage/profile/truth-injection 协议下冻结，我们没有把 CycSim "
        "事后塞进旧面板做非等价头对头，而是保留设计层比较，并把匹配条件下的"
        "头对头实验列为后续工作。",
        "",
        "线程问题在现有代码中无法复现。当前生成器提供 process-level workers、"
        "参数校验、单元测试和 1/2/4-worker 扩展性 benchmark。如果仍有具体"
        "报错，需要保留原始命令、操作系统和 release tag，再按那个版本定位。",
        "",
        "## 当前边界",
        "",
        "六物种预注册 gate 仍未通过，因此全文不写成所有物种普遍优于现有"
        "工具。当前支持的表述是：GC 条件化在具有明显局部 GC 结构的基因组"
        "中产生可测的方向性效应，受控 1-bin/5-bin 对比和拟南芥独立复现用于"
        "支撑机制层面，人类 chr21 30x 的 R3/R5 排名用于支撑真实数据锚定层"
        "的一致性。",
        "",
    ]
    args.output.write_text("\n".join(lines), encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
