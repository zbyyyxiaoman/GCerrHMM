# BMC Bioinformatics Collection: GC-Aware Error Modeling via HMM

## 技术路径与实验设计方案 v2.0 (大基因组版)

> 目标期刊: BMC Bioinformatics (JCR Q1 / 新锐二区)  
> Collection: Simulation Methods and Benchmarking in Bioinformatics  
> 投稿截止: 2026年9月15日  
> 更新: 2026-08-12 - 根据最新运行日志刷新进度与剩余目标（见第十节）

---

## 一、物种选择 (更新版)

### 1.1 六级梯度设计 (从简单到复杂)

| 层级 | 物种 | 基因组大小 | GC% | 分类 | 角色 | 优先级 |
|------|------|-----------|-----|------|------|--------|
| **L1** | *E. coli* K-12 | 4.6 Mb | 50.8% | 细菌 | 流程验证/基线 | P0 |
| **L2** | *S. cerevisiae* S288C | 12.1 Mb | 38.2% | 真菌 | 真核模式生物 | P0 |
| **L3** | *A. thaliana* | 135 Mb | 36.0% | 植物 | 中等复杂度 | P1 |
| **L4** | *D. melanogaster* | 180 Mb | 42.0% | 昆虫 | 动物模式生物 | P1 |
| **L5** | *M. musculus* GRCm39 | 2.7 Gb | 42.0% | 哺乳动物 | 大基因组 | P2 |
| **L6** | *H. sapiens* GRCh38 | 3.1 Gb | 40.9% | 人类 | 超大基因组/黄金标准 | P2 |

### 1.2 大基因组处理策略

对于人类(3.1Gb)和小鼠(2.7Gb)，提供两种运行模式：

| 模式 | 处理方式 | 适用场景 | 内存需求 | 时间 |
|------|---------|---------|---------|------|
| **染色体级** | 提取单条染色体(chr21/chr19) | 快速验证/调试 | 32-64GB | 4-8h |
| **分块全基因组** | 100Mb窗口并行 | 正式实验 | 128-256GB | 2-3d |
| **完整全基因组** | 一次性处理 | 最终验证 | 512GB+ | 3-5d |

### 1.3 推荐递进路线

```
Day 1-2:   L1-L2 (E.coli + Yeast)     → 验证流程完整性
Day 3-5:   L3-L4 (Arabidopsis + Fly)   → 验证真核生物/跨域泛化
Day 6-8:   L5-L6 chr级 (chr21/chr19)   → 验证大基因组处理能力
Day 9-14:  L5-L6 全基因组              → 生产级完整数据
```

---

## 二、核心技术创新 (不变)

### errHMM 设计

```
States: M(Match), S(Substitution), I(Insertion), D1-D3, D4+(Deletion)
GC-Aware: P(transition | GC_bin), GC_bin = floor(GC_content × 10)
训练: 从真实BAM中提取 GC_bin → 状态转移频率 → Laplace平滑
输出: 10个 8×8 转移概率矩阵 (对应10个GC_bin)
```

### 三条对比路线

| 路线 | 方法 | 说明 |
|------|------|------|
| Route A | Sample | 直接从参考样本采样错误模式 (PBSIM3默认) |
| Route B | qsHMM | 对Quality Value建模，间接推断错误 |
| Route C | **errHMM (Ours)** | 直接对错误类型建模，显式GC-aware |

---

## 三、实验流程 (适配大基因组)

### 3.1 数据准备阶段优化

```
Phase 1 改进:
  1. 使用 ncbi-datasets CLI 下载大基因组 (比wget稳定)
  2. 支持断点续传 (wget -c)
  3. 大基因组自动验证MD5
  4. 染色体提取工具 (seqkit/Custom Python)
```

### 3.2 HMM训练阶段优化

```
Phase 2 改进 (大基因组适配):
  原方案: 加载全基因组到内存 → GC计算
  新方案: 
    a) 流式处理BAM (pysam迭代)
    b) 分染色体训练 → 合并统计
    c) 增量式GC_bin计数 (不需要全加载)
  
  内存优化:
    - 不保存参考序列字典，流式查询
    - GC窗口预计算为numpy数组 (memory-mapped)
    - 转移矩阵使用稀疏表示
```

### 3.3 模拟生成阶段优化

```
Phase 3 改进 (大基因组适配):
  - 分块生成: 100Mb区块并行
  - 输出流式压缩: pigz实时压缩
  - 进度检查点: 每1000条read保存状态
  - 支持从检查点恢复 (防中断)
```

### 3.4 下游任务优化

```
Phase 5 改进 (大基因组适配):
  Mapping:
    - minimap2分索引 (每1Gb一个索引块)
    - 流式输出SAM直接转BAM
  
  Assembly (关键!):
    - <50Mb: 本地Titan运行
    - 50-500Mb: Rhea节点运行
    - >500Mb (人/鼠全基因组): Rhea节点 + 外存模式
    - hifiasm参数调优: -f0 (跳过纠错, 模拟数据已高质量)
  
  Variant Calling:
    - DeepVariant分窗口并行
    - 结果合并 (bcftools concat)
```

---

## 四、三节点协同策略 (更新)

### 节点角色

| 节点 | 角色 | 处理范围 | 大基因组特殊任务 |
|------|------|---------|----------------|
| **Titan** | CPU主控 | ≤200Mb基因组全流程 | HMM训练、模拟生成、比对 |
| **Rhea** | 内存/拼接专精 | >200Mb基因组拼接 | hifiasm大基因组拼接 |
| **Moon** | I/O/下载 | 所有数据下载 | ncbi数据集下载、格式转换 |

### 大基因组专用调度

```bash
# 人类全基因组拼接 (必须在内存节点上)
ssh <compute-node> "
    cd \"$PROJECT_DIR\" &&
    ulimit -v 41943040 &&          # 限制虚拟内存40GB
    hifiasm -o asm_human -t 16 \\
        -f0 \\
        -z20 \\
        --hg-size 3.1g \\
        data/simulated/route_C_errhmm/Hsapiens/*.fastq.gz
"

# 小鼠全基因组拼接
ssh <compute-node> "
    cd \"$PROJECT_DIR\" &&
    hifiasm -o asm_mouse -t 16 \\
        -f0 \\
        --hg-size 2.7g \\
        data/simulated/route_C_errhmm/Mmusculus/*.fastq.gz
"
```

---

## 五、资源估算

### 5.1 磁盘空间

| 阶段 | quick_test | chr_level | full_genome |
|------|-----------|-----------|-------------|
| 参考基因组 | 0.5 GB | 3.5 GB | 6.5 GB |
| 真实测序数据 | 5 GB | 40 GB | 200 GB |
| HMM模型 | 0.1 GB | 0.5 GB | 1 GB |
| 模拟数据 | 10 GB | 80 GB | 400 GB |
| 比对结果 | 5 GB | 40 GB | 200 GB |
| 拼接结果 | 2 GB | 15 GB | 100 GB |
| 结果/图表 | 1 GB | 5 GB | 10 GB |
| **总计** | **~24 GB** | **~184 GB** | **~918 GB** |

### 5.2 内存需求

| 任务 | E.coli | Human chr21 | Human Full |
|------|--------|-------------|------------|
| HMM训练 | 4 GB | 16 GB | 32 GB |
| 模拟生成 | 4 GB | 8 GB | 16 GB |
| minimap2比对 | 4 GB | 16 GB | 64 GB |
| hifiasm拼接 | 8 GB | 64 GB | **512 GB** |
| 变异检测 | 4 GB | 16 GB | 32 GB |

### 5.3 时间估算 (8线程)

| 任务 | E.coli | Human chr21 | Human Full |
|------|--------|-------------|------------|
| 数据下载 | 10 min | 2 h | 8 h |
| HMM训练 | 30 min | 3 h | 12 h |
| 模拟生成 | 1 h | 4 h | 16 h |
| 比对 | 30 min | 3 h | 12 h |
| 拼接 | 2 h | 12 h | 48 h |
| 变异检测 | 1 h | 6 h | 24 h |
| **单物种总计** | **~5 h** | **~30 h** | **~120 h** |

---

## 六、运行模式选择指南

```
你是第一次跑这个流程？
  └─ 是 → quick_test (2-4小时确认流程通畅)
  └─ 否 → 你有真实的人类/小鼠测序数据？
           └─ 否 → chr_level (用模拟数据演示)
           └─ 是 → 磁盘空间 > 1TB ?
                    └─ 是 → full_genome (完整论文数据)
                    └─ 否 → chr_level (优先完成核心结果)
```

### 快速启动

```bash
# 1. 首次配置
ssh <login-node>
bash scripts/01_setup_env.sh
source ~/.bashrc

# 2. 快速测试 (推荐第一次)
screen -S bmc_test
bash scripts/03_master_pipeline_v2.sh --mode quick_test

# 3. 确认OK后，跑正式实验
bash scripts/03_master_pipeline_v2.sh --mode chr_level

# 4. 最终全基因组
bash scripts/03_master_pipeline_v2.sh --mode full_genome
```

---

## 七、风险评估与应对

| 风险 | 概率 | 影响 | 应对策略 |
|------|------|------|---------|
| 人类基因组下载失败/太慢 | 中 | 延迟 | 使用ncbi-datasets CLI或国内镜像 |
| Rhea节点内存不足(拼接) | 高 | 任务失败 | 换用flye(内存友好)或分块拼接 |
| 磁盘空间不足 | 中 | 中断 | 监控df -h，及时清理中间BAM |
| 模拟数据质量不达标 | 低 | 需重跑 | 先跑quick_test验证 |
| PBSIM3集成失败 | 低 | 使用自研模拟器 | generate_simulated.py已可独立运行 |

---

## 八、与v1版本的主要变化

| 维度 | v1 (小基因组) | v2 (大基因组) |
|------|--------------|---------------|
| 物种 | 6个小基因组 | 6个含人类/小鼠 |
| 最大基因组 | 47 Mb (human chr21) | 3.1 Gb (human full) |
| 运行模式 | 单一模式 | 三级模式 (quick/chr/full) |
| 拼接节点 | Titan | Rhea (大内存) |
| 磁盘需求 | ~50 GB | ~200-900 GB |
| 染色体提取 | 无 | 支持chr21/chr19等 |
| 断点续传 | 无 | 支持 |
| 渐进式路线 | 无 | 4级递进验证 |

---

## 九、里程碑检查点

| 检查点 | 目标 | 验证方式 | 实际状态 (2026-08-12) |
|--------|------|---------|----------------------|
| CP1 (Day 1) | E.coli流程跑通 | `ls results/tables/level1_Ecoli_*.json` | ✅ 完成 |
| CP2 (Day 2) | 6物种参考基因组就绪 | `cat PROJECT_INDEX.md` | ✅ 完成 |
| CP3 (Day 4) | HMM训练完成 | `ls data/trained_models/*.json \| wc -l` | ⚠️ 人/鼠模型为旧错误数据训练，待重训 |
| CP4 (Day 6) | 90组模拟数据生成 | `find data/simulated -name "*.fastq*" \| wc -l` | ⚠️ quick_test 36组已齐；人/鼠需随重训再生成 |
| CP5 (Day 8) | Level-1评估完成 | `ls results/tables/level1_*.json \| wc -l` | ⚠️ 30/36（缺Dmelanogaster，无真实ONT fastq） |
| CP6 (Day 12) | 下游任务完成 | `ls results/tables/mapping_*.json \| wc -l` | ⚠️ 36组表已齐，但人/鼠为旧数据产物 |
| CP7 (Day 14) | 图表+报告生成 | `ls results/figures/*.png` | ⚠️ 旧版(08-10)，缺Dmelanogaster与fig6/fig7 |
| CP8 (Day 16) | 论文初稿 | 人工撰写 | ⬜ 未开始 |

---

## 十、当前进度与剩余目标 (2026-08-12，依据运行日志)

### 10.1 进度实况

**已完成：**
- Ecoli / Scerevisiae / Athaliana / Hsapiens_chr21 / Mmusculus_chr19 五物种 quick_test 全流程（level1/mapping/variant/assembly 各 30 组，08-10）。
- Dmelanogaster phase2-5 于 2026-08-12 14:47 跑完：模拟 6/6，mapping/variant/assembly 各 6 组（Rhea 提交失败，已回退 Titan 本地拼接）。
- 正确的人/鼠 ONT 数据已下载并校验：`Hsapiens_ont.fastq.gz` = ERR13491992、`Mmusculus_ont.fastq.gz` = SRR14685232（08-10）。

**未完成 / 存疑：**
- 人/鼠训练链路仍是旧数据产物：现存 `*_train_aligned.bam`、`Hsapiens_chr21/Mmusculus_chr19` 模型与模拟结果均为 08-09 前由错误 reads（Pseudomonas / L. reuteri）生成；08-12 10:57/11:11 的重新比对只有 START 无 DONE 记录，bam 时间戳未更新。
- `Hsapiens_chr21_ont.fastq.gz` / `Mmusculus_chr19_ont.fastq.gz` 仅 0.7/1.7 MB，是从错误数据中提取的，需从正确 fastq 重新提取。
- Dmelanogaster Level-1 被跳过：`Dmelanogaster_ont.fastq.gz` 缺失（`DRR664371.sra` 已下载但未转换）。
- 图表与 `results/stats/FINAL_REPORT.md` 仍为 08-10 旧版，不含 Dmelanogaster；`visualization.py` 新增的 fig6/fig7 钩子尚未运行。
- Level-1 error-rate 子指标（`ENABLE_LEVEL1_ERROR_RATE=1`）尚未执行。
- `scripts/03_master_pipeline_v2.sh` 第 804 行附近存在 `;;` 语法错误，运行以 rc=2 收尾（各 phase 本身已完成）。

### 10.2 剩余目标要求（验收标准）

| 编号 | 目标 | 验收标准 |
|------|------|---------|
| G1 | 重建人/鼠训练数据链：正确 ONT → train_ref 比对 → 提取 chr21/chr19 真实 fastq → 重训模型 | 比对日志出现 DONE；`*_train_aligned.bam` 与 `data/trained_models/Hsapiens_chr21_errhmm.json`、`Mmusculus_chr19_errhmm.json` 时间戳晚于 2026-08-12 |
| G2 | 清理人/鼠旧模拟与结果，重跑 phase3-5 | `results/tables/{mapping,variant,assembly}_Hsapiens_chr21_*`、`_Mmusculus_chr19_*` 时间戳晚于 G1 完成时间 |
| G3 | 补齐 Dmelanogaster 真实数据并补跑 Level-1 | `fasterq-dump DRR664371.sra` 生成 `Dmelanogaster_ont.fastq.gz`；`ls results/tables/level1_Dmelanogaster_*` = 6 个文件 |
| G4 | 受控顺序执行 Level-1 error-rate 评估 | `ENABLE_LEVEL1_ERROR_RATE=1` 全物种跑通且无 OOM |
| G5 | 重跑 Phase 6-7（消融/可视化/最终报告） | `results/figures/*.png` ≥ 7（含 fig6/fig7）；`FINAL_REPORT.md` 覆盖 6 物种且时间戳刷新 |
| G6 | 修复 `03_master_pipeline_v2.sh` 语法错误 | 脚本 `bash -n` 通过，完整运行 rc=0 |
| G7 | （可选）按磁盘/内存余量决定是否升级到 chr_level / full_genome 模式 | 磁盘 > 1 TB 且 G1-G5 全部通过后评估 |
| G8 | 修复 variant 评估：Phase 3 前加"参考基因组变异化"（如 0.1% SNP + 少量 indel），输出 truth VCF，流水线传 `--truth` 激活已写好的 bcftools norm + isec F1 代码（downstream_tasks.py:155-190）；同时核查 Ecoli route_C indel 率异常低（0.108%/base，仅为 Dmelanogaster 的 1/13，疑似模型训练数据源有误） | variant 结果含 precision/recall/F1；模拟数据 SNP F1 ≥ 0.8；truth VCF 数量/类型/生成参数留档备 Methods 使用 |

**执行顺序：** G6（先修脚本）→ G1 → G2 → G3 → G4 → G5 →（G7 评估）。G8 涉及 Phase 3 改动，应在 G1/G2 重跑之前并入，避免人/鼠重复跑两遍。

> G8 背景：variant "全 0" 已定性为设计缺陷（模拟 reads 与参考同源，无变异可检）。此前的临时修法（cigar_metadata 计数）只是绕过问题，不构成下游任务评估，论文中不可作为 variant calling 结果报告。备选：若不做预植变异，则应把 variant 从下游任务撤下，cigar 计数降级为 Level-1 indel 率一致性检查。

---

*文档版本: v2.1 | 更新日期: 2026-08-12（依据 logs/master_v2.log、logs/run_species_Dmelanogaster.log 等最新日志刷新）*
