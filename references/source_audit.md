# OE资料审查记录

## 完整性与阅读层次

已对原目录全部 271 个文件建立路径、大小、SHA-256、格式、作用、数据流与复用建议索引，约 5.03 GB。`audit/file_index.tsv` 是逐文件入口；`audit/code_inventory.json` 记录 55 个代码/说明文件的函数、读写接口和静态问题位置；`audit/code_review.tsv` 给出逐源码的实质用途与修订建议；`audit/duplicates.json` 标记 56 组字节完全相同的文件。

审查重点为源代码、56 页课程 PDF 文本、GEO 元数据及表格/矩阵结构。二进制模型权重、安装程序、RData 和大型矩阵以索引/格式/关联代码说明作用；没有为了“阅读”而运行安装器、执行打包 Perl 或载入 DrugReflector 权重。课程数据的实际运行验证另见 validation.md。此记录不声称对全部原始测量值或蛋白结构逐项作了科学验证。

## 模块取舍

| 模块 | 原项目用途 | 本 Skill 处理 |
|---|---|---|
| 01、03–09 | Alzheimer bulk GEO 注释、标准化、ComBat、limma | 总结思路，独立于单细胞扰动 |
| 10–14 | GSE163577 单细胞/单核分析、注释、CellChat | 保留 QC/定位思想，重做输入和供体契约 |
| 15–20 | eQTL/GWAS MR、交集、森林图、列线图 | 仅作扩展；指定靶点不以 MR 筛选为前置 |
| 21–26 | bulk 高低表达 GSEA/GSVA、免疫浸润与外部验证 | 观察性补充，不代表干预后状态 |
| 27–31 | 基因分布、GSN KO/OE、富集 | 本次核心重构 |
| 32–35 | Python、v-score、DrugReflector、排名图 | 独立药物候选模块 |
| 36–37 | GSN PDB、对接截图 | 独立结构/对接展示，不证明调控关系 |

## 会改变结论的修正

1. `12.Seurat/geo12.Seurat.R` 按目录前缀推断疾病组，并设 `pbmc$Sample=pbmc$Type`（约190行）。Sample 实际成了疾病状态，不能用于患者重复。本版强制显式 donor_id，状态与个体分开。
2. 原 `28`、`30` 只筛 `Type==Disease`，没有按细胞类型建网；容易由细胞组成差异制造依赖关系。本版同类同状态建网，并核查供体。
3. 原特征集只有 target+5000HVG，PPARA 可能被丢弃。本版主动核验读出基因，保留合格读出；缺失/低检出与不显著分开输出。
4. 原代码选 HVG 后再让引擎计算 QC，线粒体比例和文库大小可能因丢基因而改变。本版先在完整计数上 QC，再选模型基因，关闭引擎重复 QC。
5. `as.data.frame(countMat[...])` 过早稠密化；新版本 Knk 也明确拒绝 data.frame 输入。本版稀疏导入、内存预估，仅在算法必要阶段使用稠密邻接。
6. 原 OE 使用非导出 `scTenifoldKnk:::` 接口，版本未锁。当前官方 API 已变化：Knk 2.0.0 QC 参数改名、默认加 CPM、默认计算 direction。新流程锁提交/函数指纹、关闭 direction、复用同一 WT。
7. 原图中的 FC/Z 为网络统计，不是表达 FC；所谓“火山图”不能用左右侧代表上/下调。新图明确为网络距离/排名与模型 FDR。
8. `29/31.Enrich` 的 `p.adjustFilter=1` 实际按 nominal P 筛选；未指定实际建网 universe；KEGG 转换引用 `rt$genes/rt$entrezID`，但 sigDiff 只有 gene/distance/Z/FC/P。本版保存显式映射、背景和 BH 结果，空集也有状态。
9. `21.GSEA` 按 raw P<0.05 取结果；`22.GSVA` 对多个通路 t 检验未校正。两者是 bulk 靶点高低表达分组，不是 KO/OE 生成的表达响应。主流程不将它们混入扰动方向证据。
10. `16.MR` 先按结局 P>0.05 筛 SNP（约61行），可能产生选择偏差；`17` 将五方法一致及多效性 P>0.05 作筛选，不能由不显著证明没有多效性。MR 未作为本次核心实现。
11. `06.sva` 第二轮循环仍用上一轮 `file` 作排除判断（约39–40行），存在错误跳过/包含输出文件风险；ComBat 未给疾病设计矩阵，需审查病情与批次混杂。旧脚本未直接运行。
12. `23.CIBERSORT` 在约133行含 `if(substr(Sys.Date(),1,4)>2026){next}`，函数内该位置不在循环，后续年份可能报错。这是课程代码的时间条件；本版不复制该模块。

## 药物模块的真实含义

`33.vscore` 以 bulk Treat-Control 均值差除以两组方差和的平方根，生成一条疾病签名（1×gene h5ad），不是细胞级表达数据。`34` 使用三个神经网络权重对化合物分类 logits 取平均并 softmax、排序；prob 不是临床有效概率，也不是虚拟 KO/OE 后 PPARA 的响应概率。需要另查签名方向、训练任务和适用域才能声称“逆转疾病”。

`drug_reflector.py` 的字符串清洗不是真正的 Ensembl→HGNC 映射；`models.py` 使用较新的 Python 类型标注，与 pyproject 的 >=3.8 声明不完全一致；顶层与内部包的 CLI/import 结构需要单独包装检查。`predict.py` 的 background-samples 参数未用于计算 P 值，`cli.py` 的 device 参数没有接入实际预测；`35.vol` 画的是 logit–softmax 概率，名称虽叫 vol，但不是差异统计火山图。signature_refinement 的相关性融合也不是因果预测。

## 副本与版本

`05.normalize/geo05.normalize.R` 与 GSE5281 副本完全相同；其他队列副本主要改变路径、GEO ID 等参数，详细差异另见 `audit/variant_diffs.txt`。DrugReflector 的 build/lib 下多个文件与工作源码完全相同，使用工作源码作主副本，避免把打包输出当第二套算法。29/31 富集脚本的少量区别主要是目录、注释/排版，算法问题相同。

源材料保持原样；新实现置于单独 Skill 目录。源码审查中发现的问题不表示原课程的每个结果必然错误，而是不能把教学脚本未经核验直接迁移到肝纤维化问题。
