# 运行环境与恢复

本机验证使用 R 4.5.1，Seurat 5.4.0，SeuratObject 5.3.0，Matrix 1.7.4，scTenifoldKnk 2.0.0，scTenifoldNet 1.4.3。Python 转换使用 Anaconda Python 3.13、h5py 3.15.1、scipy 1.16.3、pandas 2.3.3；详细版本记录在 `dependencies.lock.json` 和验证运行的 sessionInfo。

资源检测：约 31.84 GB 总内存、检测时约 11.77 GB 可用、16 逻辑 CPU；可用内存随其他程序变化。核心网络使用 CPU，GPU 不是必需。

`references/environment.local.json` 保存本机 Rscript/Python 路径；`library_paths.local.txt` 保存项目独立 R 库及已有包库。个人技能目录只装源码和说明；依赖仍引用当前工作区 `.runtime/r-library`。**搬移或删除工作区后需更新路径或重建独立库。** 本机配置文件由安装者创建，不进入 Git 仓库；可参考 `environment.example.json`。

## 新机器或库损坏

1. 安装适配平台的 R；按 `dependencies.lock.json` 恢复依赖版本。核心依赖为 Matrix、MASS、RSpectra、cli、ps、RhpcBLASctl、enrichR、igraph、reshape2；接口另需 Seurat、SeuratObject、yaml、jsonlite、digest。GO 可选 clusterProfiler 与物种 OrgDb；测试需 testthat。
2. 使用新的项目库，不在本流程中升级其他项目的全局包。将已有包库写入 `library_paths.local.txt` 或环境变量 `SCVP_R_LIBS`（分号分隔）。缺包时按报错补齐，再安装固定源码。
3. `Rscript scripts/install_engine.R --lib <项目库>` 会先核对所附上游文件 SHA-256，再安装固定版本。这些快照只含构建所需 R 源码和命名空间，不含上游完整帮助页。原作者与 GPL 授权见 DESCRIPTION / `LICENSE`。
4. 在新的 R 进程运行 `Rscript scripts/run_tests.R`，以及 `python tests/test_h5ad.py`。改变版本后须重新核验实际函数行为，不能直接修改 engine-lock 让检查通过。
5. 运行最小课程演示或经授权的小数据测试，再处理正式全量数据。

Windows R 的中文路径由 bootstrap 设置 UTF-8 locale；Python Matrix Market 使用文件句柄绕过当前 scipy 路径兼容问题。上游源文件保持原样。初始 R 启动可能因宿主 `C.UTF-8` locale 出现警告，bootstrap 后解析和实际运行已验证。

## 验证命令

```powershell
Rscript --vanilla scripts/run_tests.R
python -X utf8 tests/test_h5ad.py
```

准备课程数据示例（源目录在技能包之外）：

```text
python scripts/prepare_course_demo.py --source <OE资料绝对路径> --output <验证目录>/course_data --config <验证目录>/course_demo.yaml
Rscript scripts/run_perturbation.R --config <验证目录>/course_demo.yaml
```

准备脚本从本机原材料抽取 4 个独立 AD 供体各至多 240 个课程注释的内皮细胞，保留所有基因后再 QC，使用 GSN 扰动，以 DLC1 和 PPARA 作读出。它会生成小型 10X 文件，不修改原材料。该演示用于流程检查，不能提供肝纤维化结论。
