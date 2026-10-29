# 方法定义、来源与证据边界

## 固定源码

- scTenifoldKnk **2.0.0**，提交 `e003c72d3a153ec07512271c1ae0fbeeadc319ca`：[源码](https://github.com/cailab-tamu/scTenifoldKnk/tree/e003c72d3a153ec07512271c1ae0fbeeadc319ca)。
- scTenifoldNet **1.4.3**，提交 `3103a53216c347a3d4c9660d63b03c2ca99b47a1`：[源码](https://github.com/cailab-tamu/scTenifoldNet/tree/3103a53216c347a3d4c9660d63b03c2ca99b47a1)。
- 本地保存对应 R 源码、DESCRIPTION、NAMESPACE；`upstream/manifest.json` 保存 URL、提交和 SHA-256。核心函数的安装后指纹在 `engine-lock.json`。
- Net 原始方法 DOI：[10.1016/j.patter.2020.100139](https://doi.org/10.1016/j.patter.2020.100139)。解释以本次固定版本源码为准，不能假定课程代码对应当前 API。

## 统计量

流形对齐为每个基因得到 WT 和扰动网络中的坐标。两套坐标的欧氏距离是 `distance`；变换后的标准化距离形成 `Z`。固定版本 `dRegulation.R` 使用：

`E = mean(distance[非目标基因]^2)`；`FC = distance^2 / E`。

随后使用自由度为 1 的卡方分布上尾得到理论模型 P 值，并用 BH/FDR 校正。该 FC 不含表达方向。它不是表达量比值、不是 log2FC。Z 低于零也不说明表达下调。P 值刻画此模型内相对网络扰动，不能取代以供体/实验重复为单位的干预效应检验。

完全相同网络直接输出 `distance=0、p=1、p.adj=1、FC/Z/rank=NA、unchanged_network`，避免除零或随机特征向量把纯数值差异转成显著结果。普通变化使用官方统计函数；数值异常显式失败。

## 为什么 OE 不能按真实表达倍数解读

原 `refer.scTenifoldOE.R` 的第一方法先建立 WT，经 strictDirection/对角线清零/转置，再执行 `OE[gOE,] <- WT[gOE,] * oe_factor`。这改变的是模型边权，没有生成实验后的细胞表达矩阵。正式交付名称为 `OE_edge_scaling`。

第二方法对目标 counts 的每个细胞加一后重新建两个网络。固定 pcNet 源码逐基因去均值、除标准差，因此直接对它输入的数据加常数或乘正数，会被标准化消除；自动化数值测试已验证加一不改变直接 pcNet 输出。完整 scTenifoldNet 在此前还对两组分别做 CPM，添加目标 counts 会改变全细胞总量并使其他基因也改变归一化值。由此产生的网络差异可能是组成效应，加上两个网络的随机抽样差异，不能视为特异性过表达响应。故未开放 `scTenifoldOE2` 为可运行模式。

原 OE 第一方法的 `makeNetworks` 直接接收 counts，而当前官方 KO 包会先 CPM；原两条路线还可能在 QC/归一化上不匹配。本版通过复用官方 KO 返回的 WT，使 KO 和 OE 基线一致；结果不承诺重现课程的原数字。

## 2.0.0 的方向预测为何关闭

新函数 `knockoutDirection` 使用 WT log1p CPM 基因相关矩阵的热核传播估计方向。其源码文档明确指出：在部分细胞系中的表现高于随机，但常反映多个 knockdown 共享的主表达程序，未必靶点特异，且可能被测序深度驱动，导致几乎所有基因预测为下降。

这属于基于共表达假设的方向模型，不是已为 FABP1–PPARA 人肝纤维化场景验证的表达响应。当前配置固定关闭 `dr_direction`；`expression_prediction.enabled=true` 会明确拒绝运行，要求先增加适配数据和验证过的实现，而不会自动给出看似精确的表达量。

## 其他限制

- FABP1 和 PPARA 的观察相关性可由肝细胞分区、脂质状态、病因、处理及供体因素造成。相关或网络连接不证明直接因果调控，也不证明调控方向。
- snRNA 与 scRNA 的检出和 QC 机制不同；固定阈值只是配置起点。某基因低检出时，增加细胞可能有帮助，插补不能当成真实测量补齐证据。
- 选定基因集影响 CPM、网络和 P 值参考分布。跨网络比较必须检查共同基因和基线差异。
- GO 可使用固定 org.Hs.eg.db 注释或显式 GMT；KEGG 要求快照。课程 GMT 用于技术演示，其当前生物学时效未验证。
- PPARA 签名评分采用供体内聚合后、同细胞类型跨供体/状态的基因 z 分数加权均值，排除 FABP1/PPARA 自身，记录覆盖率。它是观察性目标程序代理，不能称作直接测得的 TF 活性或模拟后活性。
