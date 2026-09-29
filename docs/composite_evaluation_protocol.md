# Composite 评估口径与公平性说明

## 当前 composite

当前 Level-1 composite 由以下维度加权：

| Metric | Weight |
|---|---:|
| Read length similarity | 0.15 |
| Error-rate similarity | 0.25 |
| QV histogram Pearson | 0.20 |
| GC KS similarity | 0.25 |
| k-mer frequency similarity | 0.15 |

当 `ENABLE_LEVEL1_ERROR_RATE` 未开启时，error-rate 项不参与计算，
其余权重按实际可用维度重新归一化。当前 cross-tool 表就是这种情况，
因此表内只出现 read length、QV、GC 和 k-mer 四项。

## 是否偏袒 PBSim3

会产生口径优势，但不能简单称为作弊：

1. PBSim3 使用的是 sample-based 方法，训练输入与评测数据来自同一
   真实读集合。它天然更容易匹配 read-length 和 QV 分布。
2. errHMM 分支虽然也已传入 profile JSON，但它的核心贡献是
   GC 条件错误状态建模，而不是直接拟合每个分布。
3. 当前 composite 中含有 read length、QV 和 GC 三个分布型指标，
   又没有启用 error-rate 项；因此它更接近“分布相似度”评分，
   而不是“生物真实性”或“下游任务价值”评分。
4. PBSim3 在 Ecoli 的 QV 子分数达到 `100.0`，与其采样拟合机制一致。

结论：composite 可用于描述同一个输入、同一个深度下的分布相似度，
但不应被解释成 PBSim3 在所有科学维度上全面优于 errHMM。

## 推荐报告方式

1. 主表保留 composite，但必须同时给出各子指标，不能只报总分。
2. 将 PBSim3 明确标注为 `sample/profile-fitted baseline`。
3. 在论文中把 composite 命名为 Level-1 distributional similarity，
   不要写成 overall simulator accuracy。
4. 单独报告 downstream mapping / variant / assembly 结果，不把其
   混入当前 composite 后再宣传统一胜负。
5. 对 error-rate 项做敏感性分析：至少补充一次启用
   `ENABLE_LEVEL1_ERROR_RATE=1` 的结果；若成本过高，则在附录明确
   说明当前 composite 未覆盖该维度。
6. k-mer 必须使用本次修复后的 matched-sampling 口径，并同时报告
   Jaccard 与 frequency correlation，避免一个指标掩盖真实的集合重合。

## 推荐论文表述

可以写：

> Under matched training data and depth, errHMM and PBSim3 exhibit
> complementary strengths: PBSim3 is strongest on distributional
> similarity, whereas errHMM is competitive on sequence-level k-mer
> fidelity and GC-aware metrics.

不应写：

> errHMM outperforms all existing simulators on the composite score.
