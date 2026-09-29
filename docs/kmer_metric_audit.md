# NanoSim Scerevisiae k-mer 指标审计

## 结论

NanoSim Scerevisiae 原始结果中的 `kmer=0.0` 不是“完全没有可用 k-mer”，
而是评估器的采样口径不一致导致的假阴性。

原始评分路径中：

1. 真实数据的 k-mer 缓存只保留 top 1,000 个 k-mer。
2. 模拟数据在比较时仍使用完整 Counter。
3. 两侧的集合规模、频率归一化口径不同，交集很容易萎缩，
   频率 Pearson 也可能被少数样本拉到接近 0 或负值。

因此原始 `kmer_freq_corr=-0.00219` 被 `max(0, ...)` 截断成 `0`。

## 直接复核

用同一套 21-mer、固定 reservoir 抽样和归一化频率，直接复核
Scerevisiae 真实 FASTQ 与各模拟 FASTQ：

| Tool | Reads | Jaccard | Frequency Pearson | Common unique 21-mers |
|---|---:|---:|---:|---:|
| NanoSim | 250 | 0.02110 | 0.15850 | 31,213 |
| NanoSim | 1,000 | 0.04442 | 0.31850 | 241,518 |
| NanoSim | 5,000 | 0.18512 | 0.32773 | 3,594,797 |
| errHMM | 1,000 | 0.03999 | 0.42532 | 203,675 |
| PBSim3 | 1,000 | 0.00960 | 0.06396 | 29,175 |
| badread | 1,000 | 0.02421 | 0.06741 | 126,015 |

NanoSim 在 1,000 reads 同口径下已经有明显的公共 k-mer 和正相关，
所以“k-mer=0”不是模拟器真实缺失。

## 修复

`src/evaluate_level1.py` 已改为：

- 真实和模拟两侧都使用相同大小的 top-N k-mer 样本；
- 两侧统一按各自总 k-mer 数归一化；
- 频率相关性和 Jaccard 在相同的 top-N 集合口径上计算；
- 缓存版本提升为 v4，避免继续读取旧 top-1,000 缓存。

修复后 Scerevisiae NanoSim：

`kmer 0.00 -> 29.78`，`composite 70.98 -> 76.85`。

所有 12 个 cross-tool 结果均已用新口径重算。
