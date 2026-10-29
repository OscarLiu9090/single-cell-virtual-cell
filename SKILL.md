---
name: single-cell-virtual-perturbation
description: 基于原始单细胞计数和供体信息，运行 scTenifoldKnk 虚拟敲除及共享 WT 网络的边权过表达模拟，评估指定下游基因的网络扰动、稳定性和供体表达关联。适用于 FABP1–PPARA 等基因对；不把网络评分解释为表达升降。
---

# 单细胞虚拟扰动

将一个扰动基因与一个或多个关注基因，在指定细胞类型和状态中做可复核分析。默认使用已锁定并验证的 scTenifoldKnk 2.0.0 / scTenifoldNet 1.4.3。FABP1–PPARA 是配置示例，不是预设机制。

## 首次使用

1. 阅读 [中文方法说明](references/workflow_zh.md) 和 [输入及配置](references/input_output.md)。先确认物种、原始 counts、样本/供体、疾病状态及细胞类型。能从元数据查明的信息不要重复询问。
2. 查看 [本机环境与恢复](references/environment.md) 及 `references/environment.local.json`。运行前核对资源；`engine-lock.json` 检查版本与函数指纹。不要为了通过检查而重新生成锁文件。
3. 将 `configs/fabp1_ppara_human.yaml` 复制到用户本次分析目录，填写数据、输出路径和元数据映射。不要在技能安装目录存入大型数据或分析结果。
4. 先用 `scripts/run.ps1 -Config <配置绝对路径> -Preflight`（Windows），或 `Rscript scripts/run_perturbation.R --config <配置> --preflight` 检查数据；查看输出的可评估性表，再运行同一配置去掉 preflight。
5. 汇报 `run_summary.json`、`run_status.tsv`、`readout_summary.tsv` 和 `donor_association.tsv`；引用具体输出路径。缺失、低检出、无变异、无出边、运行失败和不显著是不同状态。

若尚未选定数据，先给出本研究的筛选标准和元数据需求；只有用户需要寻找数据时再检索公开数据。没有真实数据时使用 [验证说明](references/validation.md) 中的技术演示，明确其疾病与组织来源。

## 解释规则

- KO 将 WT 的目标基因出边置零；OE 将同一 WT 的目标行乘以倍数。OE 是自定义网络敏感性模拟，2/5/10 倍边权不等于 2/5/10 倍 mRNA 或蛋白。
- `distance`、`Z`、`FC` 是无表达方向含义的网络统计量。`FC = distance² / mean(non-target distance²)`。Z 的正负不是基因上/下调；模型 P 值不是患者层面的干预效应检验。
- 固定 `dr_direction=FALSE`。2.0.0 新增的方向预测依赖 WT 共表达，官方源码也说明其特异性与文库大小限制。没有适配真实扰动数据或独立验证模型时，表达方向写“尚不能判断”。
- 供体层面的 FABP1–PPARA 相关性是观察关联。PPARA mRNA、PPARA 靶基因程序评分和网络受扰动程度分别说明；不互相替代。
- 同类细胞、同状态内建网；供体平衡抽样。种子重复和留一供体分析不增加生物学样本量。不同网络的距离和 FC 不直接作差；结合排名稳定性判断。
- 不筛选出支持预设结论的种子、倍数或细胞亚群。读出基因可评估性、无显著结果和对照扰动一并交付。

## 按任务读取

- 理解课程项目或核对旧脚本：[资料审查](references/source_audit.md)，逐文件索引 `references/audit/file_index.tsv`，逐源码审查 `references/audit/code_review.tsv`。
- 判断能否回答“表达升降”或调整 OE 模式：[方法边界](references/method_limits.md)。备用“加一”方案不开放运行。
- 更换输入类型、物种、细胞或读出：[输入及配置](references/input_output.md)。原始数据需要先完成 QC、聚类、人工核验注释；本引擎不猜测细胞身份。
- 修复或更新引擎：核对 `references/upstream/manifest.json` 中固定提交，运行 `scripts/run_tests.R` 和 `tests/test_h5ad.py`，然后重做技术演示。包更新需要显式重新验证。

## 交付

提供方法和数据说明、参数/版本、全部网络结果、关注基因摘要、供体表达表、稳定性检查、富集背景与来源、日志及失败状态。指出当前证据能回答什么，以及仍需要哪些扰动实验；不要仅提供显著基因图。
