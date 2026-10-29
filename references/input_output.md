# 输入、配置与输出

## 输入契约

矩阵为非负、有限、接近整数的原始计数。R 端统一为 **基因 × 细胞**稀疏矩阵；h5ad 的 **细胞 × 基因**在转换时转置。整数检查只能发现常见误用，不能证明数据真为原始计数，仍应核对来源。

支持：

| `input.type` | 必需配置 | 说明 |
|---|---|---|
| `seurat_rds` | `path` | RNA assay 的 counts；元数据可内置 |
| `seurat_rdata` | `path`，必要时 `object_name` | 在隔离环境 load，多个 Seurat 对象时必须指定 |
| `10x` | `samples_file` | 原始 mtx/features/barcodes，支持 gzip 和未压缩；不同样本基因集取并集 |
| `h5ad` | `path`，`counts_layer` 默认 counts | 仅明确指定的 `layers/counts`；不猜测 X/raw.X；需要 Python h5py/scipy/pandas |
| `matrix_rds` | `path`、`metadata` | 带行列名的 R 矩阵，便于与现有工具衔接 |

Seurat v5 多 counts 层在细胞互不重叠时 JoinLayers；若层内细胞重复，报错要求先选定真正原始层。使用 `assay: RNA`；拒绝把 SCT/integrated 值作为原始计数。

基因名称默认已是当前物种 symbol。Ensembl 输入需提供 `gene_mapping` TSV，列为唯一 `gene_id` 与 `symbol`；h5ad 可用 `gene_symbol_column` 读取 var 中的符号列。多个 feature 映射到同一 symbol 时合并计数并输出映射表，不用 make.unique 制造假符号；映射有缺失则停止。人用 FABP1/PPARA，小鼠用 Fabp1/Ppara，不能简单把所有物种基因转成大写。

细胞元数据 TSV 至少需要以下列；没有的信息填写 NA，不捏造。

```text
cell_id                 sample_id    donor_id   condition   cell_type
S01_AAAC...-1           S01          D01        Control     Hepatocyte
S02_TTGC...-1           S02          D02        Fibrosis    Hepatocyte
```

`sample_id` 是取样/文库，`donor_id` 是人或动物生物个体。同一供体多次取样不能当多个独立供体。缺 donor 会跳过建网；缺 condition/cell_type 停止并要求完成元数据或注释。缺 sample_id 可记录为 NA，但要核实原始映射。`donor_id` 不使用保留字 `none`（该值表示没有留出供体）。

外部 metadata 覆盖对象内同名列，必须覆盖所有导入细胞。通过 `metadata_columns: {cell_type: cellType, condition: Type, donor_id: donor}` 显式映射旧列名；不可把课程 `Sample=Type` 当供体。

10X 样本表：

```text
sample_id   path                  donor_id   condition
S01         ../raw/S01            D01        Control
S02         ../raw/S02            D02        Fibrosis
```

每个 barcode 加上 `sample_id_` 前缀；外部 metadata 的 cell_id 必须一致。样本表中的 path 相对于样本表目录，配置中的其他路径相对于 YAML 所在目录。

## 推荐操作

复制 `configs/fabp1_ppara_human.yaml` 到新的分析目录，设置数据和输出路径。先在读入数据的工具中看完整 QC 分布和注释，再执行：

```powershell
# 在技能目录下；也可以用绝对脚本路径。
./scripts/run.ps1 -Config 'C:/your_analysis/config.yaml' -Preflight
./scripts/run.ps1 -Config 'C:/your_analysis/config.yaml'
```

等价 R 入口：

```text
Rscript scripts/run_perturbation.R --config C:/your_analysis/config.yaml --preflight
Rscript scripts/run_perturbation.R --config C:/your_analysis/config.yaml
```

Preflight 会导入、QC、核对可评估性，并生成供体表达/关联摘要，不构建网络。每次生成独立时间戳/PID 输出目录，不覆盖旧结果。退出码：0 为流程完成（仍必须检查每组 status），1 为致命输入/流程错误，2 为没有完成网络或部分网络运行失败。

## 参数选择

- 默认人肝细胞、Control/Fibrosis、FABP1、PPARA。默认 3,000 HVG，10 子网络，每网至多 500 个细胞，3 PC、张量 K=3、流形维度=2、q=0.9。1 核默认便于复现并限制内存；不是必须串行的算法。
- 按最少供体细胞数平衡候选池，每供体上限 2,000。最少供体异常少时应检查样本，而非让全数据被其限制后仍报告强结论。
- 5 个预设种子；供体留一使用第一枚种子，至少 3 个供体才执行。正式分析建议更多供体，必要时加大重复，不把算法重复当独立样本。
- `min_target_detection=0.05` 同时用于靶点/关注基因的预检查。低检出关注基因可出现在全量建网表，但专项摘要会标记不可可靠评估，数值置为 NA。
- `memory_budget_gb=6` 是保守的网络工作空间预算。内存估算约 `8*G²*(4*nNet+20)` 字节，尚需容纳原始数据和 R；不得将该估算视为真实峰值的保证。
- `allow_exploratory_low_donors=false`；设为 true 仅允许明确标注的低供体探索，不消除数据限制。
- 原始数据的 doublet/环境 RNA 处理不由本脚本自动推断；配置可接入已审核双细胞标签。QC 阈值需按组织和测序技术选择。
- `expression_prediction.enabled=false` 是当前能力边界。true 会报 `UNVALIDATED_PREDICTOR`，不会触发不适配模型。

## GO/KEGG 与签名输入

默认 GO 使用物种 OrgDb，可配置 `go_gmt` 替代；`kegg_gmt` 必须指定有版本的快照，否则输出 `offline_KEGG_GMT_required`。GMT 每行是 `term<TAB>description<TAB>symbol...`，符号须与建网基因一致。每次输出源文件 SHA-256、实际映射背景和入选基因。

可选签名 TSV 列：`set、gene、weight、type`，每个集合内 gene 唯一；另附 JSON 至少 `source` 与 `version`，建议包含 URL、物种、检索日期和筛选方法。`type` 建议为 `regulon_proxy` 或 `pathway`。权重正负表示先验方向，不是本次模拟推导的方向。默认至少 5 个有效基因，输出总数与覆盖数。配置开关为 `signatures.enabled`。

## 输出含义

| 文件 | 用途 |
|---|---|
| `report_zh.md`、`run_summary.json` | 运行概况及解释边界 |
| `run_status.tsv` | 缺失、数据不足、失败与已完成状态 |
| `target_readout_eligibility.tsv` | 每组靶点/关注基因的检出率、变异、可评估性 |
| `readout_per_run.tsv` | 每个种子/留一/倍数的关注基因结果，未建网读出保留 NA 状态 |
| `readout_summary.tsv` | 排名范围、显著比例、有效状态；方向固定 not_estimated |
| `stability_all_genes.tsv`、`seed_rank_agreement.tsv` | 全基因稳定性和种子间排名一致性 |
| `donor_expression.tsv`、`donor_association.tsv` | 按供体聚合的原始计数/CPM与状态内表达关联 |
| `donor_signature_scores.tsv`（可选） | 观察性靶基因程序评分 |
| 组/holdout/seed 下 `*_all.tsv`、`*_significant.tsv` | 全量模型结果和校正后显著结果 |
| `WT.rds`、`sampled_cells.tsv`、`network_genes.tsv` | 可复核的基线、候选细胞、基因与匹配对照 |
| `enrichment/`、`enrichment_status.tsv` | GO/KEGG 全量与显著结果、背景、映射及来源 |
| `input_summary.json`、`config.yaml`、`sessionInfo.txt`、`installed_packages.tsv` | 输入计数/元数据指纹、参数和版本 |
| `run.log`、`warnings.txt` | 运行诊断；无出边、数值与可评估性限制需报告 |

不得只读取显著表就把缺失的 PPARA 解释为“不受影响”。先查 eligibility/status，再读全量与稳定性摘要。
