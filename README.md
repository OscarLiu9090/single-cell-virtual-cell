# Single-cell virtual perturbation

可调用的 Codex Skill：单细胞虚拟敲除、共享 WT 的网络边权过表达模拟、指定下游基因的扰动排名及供体表达关联。附人肝纤维化 FABP1–PPARA 配置。

**网络 FC/Z 不代表 mRNA 表达倍数或方向。** 本版不会从未经验证的共表达模型输出“PPARA 必然升高/降低”。

## 内容

- [Skill 入口](SKILL.md)
- [中文方法说明](references/workflow_zh.md)
- [资料代码审查](references/source_audit.md)与[逐文件索引](references/audit/file_index.tsv)
- [输入、配置、输出](references/input_output.md)
- [方法边界与源码版本](references/method_limits.md)
- [环境恢复](references/environment.md)与[验证结果](references/validation.md)

## 使用

将仓库目录放入 `~/.codex/skills/single-cell-virtual-perturbation`，按环境文档配置 R/Python 与包库。后续可直接说：

> 使用 $single-cell-virtual-perturbation 分析我的肝纤维化肝细胞数据，评估 FABP1 KO 和边权 OE 对 PPARA 的网络影响。

或者复制配置到分析目录，填写自己的数据位置后运行：

```bash
Rscript scripts/run_perturbation.R --config /path/to/config.yaml --preflight
Rscript scripts/run_perturbation.R --config /path/to/config.yaml
```

支持 Seurat RDS/RData、10X、原始 counts 层的 h5ad，以及 R 矩阵。正式网络建模需要原始计数、供体、状态及已审核细胞类型。

核心版本：scTenifoldKnk 2.0.0、scTenifoldNet 1.4.3，固定 Git 提交、文件哈希与安装后函数指纹。课程数据完成9次缩小网络运行；人肝纤维化示例尚待实际数据。仓库不包含课程原始数据、模型权重、私人分析结果或本机环境路径文件。

## 开发验证

```bash
Rscript scripts/run_tests.R
python tests/test_h5ad.py
```

包内固定上游源码保留原作者和 GPL 授权；本项目按 GPL-3.0-or-later 提供。课程来源为用户提供的生信自学网 OE 资料，仓库提供审查与独立整理，不分发原课程数据包。
