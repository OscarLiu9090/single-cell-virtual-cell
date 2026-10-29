"""Index all original files without executing original scripts or model weights."""
import argparse
import ast
import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path

ROLES = {
    '01': ('课程说明', '项目全流程与演示图', '阅读参考'),
    '03': ('GEO 原始资料', 'GEO→04 注释', '背景数据'),
    '04': ('探针注释', '表达矩阵/平台注释→geneMatrix', 'bulk 扩展'),
    '05': ('bulk 标准化', 'geneMatrix+s1/s2→normalize', '不用于单细胞原始计数'),
    '06': ('bulk 合并去批次', '各队列 normalize→merge', 'bulk 扩展'),
    '07': ('箱线图', '合并前后矩阵→分布图', '可视化参考'),
    '08': ('PCA', '合并前后矩阵→PCA', '可视化参考'),
    '09': ('bulk 差异', '合并矩阵→limma DEG', '可选背景证据'),
    '10': ('单细胞下载元数据', 'GEO→样本信息', '技术演示来源'),
    '11': ('单细胞整理', 'GEO+临床表→12', '技术演示来源'),
    '12': ('Seurat 预处理与注释', '10X→Seurat.Rdata/细胞注释', '核心重构，供体与细胞类型需核验'),
    '13': ('CellChat', '归一化表达+细胞类型→通讯', '可选扩展，不证明细胞内 FABP1→PPARA'),
    '14': ('单细胞绘图', 'Seurat+markers→图', '可视化参考'),
    '15': ('MR 输入', 'eQTL+AD GWAS→16', '可选扩展'),
    '16': ('孟德尔随机化', '暴露/结局→MR结果', '需重新审查设计，不作为扰动前提'),
    '17': ('MR 筛选', 'MR/多效性→IVW筛选', '可选扩展'),
    '18': ('DEG 与 MR 交集', '09+17→候选基因', '已指定靶点时无需此步骤'),
    '19': ('森林图', 'MR结果+候选→图', '绘图参考'),
    '20': ('列线图模型', '表达+候选→诊断图', '无关主线'),
    '21': ('GSEA', 'bulk 靶点高低表达→排序富集', '观察关联，不代表虚拟扰动'),
    '22': ('GSVA/ssGSEA', 'bulk+基因集→通路评分', '可选观察关联'),
    '23': ('CIBERSORT', 'bulk+参考→免疫比例', '不作为单细胞前提'),
    '24': ('免疫比例图', 'CIBERSORT→图', '绘图参考'),
    '25': ('免疫相关性', 'bulk+免疫比例→相关图', '观察关联'),
    '26': ('外部 bulk 验证', 'GSE36980+候选→图', '疾病特异数据需替换'),
    '27': ('靶基因定位', 'Seurat+候选→DotPlot', '核心定位思路保留'),
    '28': ('虚拟 KO', '疾病 counts→scTenifoldKnk→网络差异', '核心重构'),
    '29': ('KO 富集', 'sigDiff→GO/KEGG', '重构背景、映射和FDR'),
    '30': ('网络边权 OE', '疾病 counts→WT→靶点行放大→差异', '核心自定义敏感性模拟'),
    '31': ('OE 富集', 'sigDiff→GO/KEGG', '同29修正'),
    '32': ('Python 安装', '安装程序/命令', '不运行旧安装器'),
    '33': ('v-score', 'bulk Treat-Control→1×gene h5ad', '药物筛选输入，非单细胞原始计数'),
    '34': ('DrugReflector', 'v-score+3个权重→化合物排名', '独立扩展，不预测 FABP1 扰动表达'),
    '35': ('药物排名绘图', 'rank/logit/prob→图', '不是差异表达火山图'),
    '36': ('对接结构准备', '蛋白结构→37', '独立扩展'),
    '37': ('对接展示', 'PDB+截图', '独立扩展'),
}


def decode(data):
    for encoding in ('utf-8-sig', 'gb18030'):
        try:
            return data.decode(encoding), encoding
        except UnicodeError:
            pass
    return data.decode('utf-8', errors='replace'), 'binary_or_mixed'


def audit(source, output):
    source, output = Path(source), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    records, code = [], []
    hashes = defaultdict(list)
    for p in sorted(source.rglob('*')):
        if not p.is_file():
            continue
        rel = p.relative_to(source).as_posix()
        h = hashlib.sha256()
        with p.open('rb') as handle:
            for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
                h.update(block)
        sha = h.hexdigest()
        role, flow, reuse = ROLES.get(rel.split('.')[0], ('其他', '', '人工判断'))
        rec = dict(path=rel, bytes=p.stat().st_size, sha256=sha, format=p.suffix.lower(), role=role, data_flow=flow, reuse=reuse)
        records.append(rec); hashes[sha].append(rel)
        if p.suffix.lower() in ('.r', '.py', '.pl', '.toml', '.md') and '__pycache__' not in rel:
            text, encoding = decode(p.read_bytes())
            entry = dict(path=rel, encoding=encoding, lines=len(text.splitlines()), functions=[], inputs=[], outputs=[], flags=[])
            if p.suffix.lower() == '.py':
                try:
                    tree = ast.parse(text)
                    entry['functions'] = [{'name': n.name, 'line': n.lineno} for n in ast.walk(tree) if isinstance(n, (ast.FunctionDef, ast.ClassDef))]
                except SyntaxError as e:
                    entry['flags'].append(f'parse_error: {e}')
            else:
                entry['functions'] = re.findall(r'([\w.]+)\s*(?:<-|=)\s*function\s*\(', text)
            patterns = {'hardcoded_workdir': r'setwd\(', 'dense_conversion': r'as\.matrix\(|as\.data\.frame\(count', 'internal_API': r':::', 'uncorrected_enrichment': r'p.adjustFilter\s*=\s*1', 'date_condition': r'Sys.Date\(', 'packed_Perl': r'package __par_pl', 'condition_as_sample': r'\$Sample\s*=\s*\w+\$Type', 'mixed_disease_cells': r'subset\s*=\s*Type\s*==\s*"Disease"'}
            for flag, pattern in patterns.items():
                lines = [i for i, line in enumerate(text.splitlines(), 1) if re.search(pattern, line)]
                if lines:
                    entry['flags'].append({'flag': flag, 'lines': lines})
            for key, pattern in [('inputs', r'(?:read\w*|load|source)\s*\([^\n]{0,180}'), ('outputs', r'(?:write\w*|save\w*|pdf)\s*\([^\n]{0,180}')]:
                entry[key] = re.findall(pattern, text, flags=re.I)[:50]
            code.append(entry)
    with (output / 'file_index.tsv').open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=records[0].keys(), delimiter='\t'); w.writeheader(); w.writerows(records)
    (output / 'code_inventory.json').write_text(json.dumps(code, ensure_ascii=False, indent=2), encoding='utf-8')
    dupes = [{'sha256': sha, 'paths': paths} for sha, paths in hashes.items() if len(paths) > 1]
    (output / 'duplicates.json').write_text(json.dumps(dupes, ensure_ascii=False, indent=2), encoding='utf-8')
    summary = {'files': len(records), 'bytes': sum(r['bytes'] for r in records), 'source_code_or_docs': len(code), 'exact_duplicate_groups': len(dupes)}
    (output / 'audit_summary.json').write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding='utf-8')
    # Keep machine-specific paths out of the portable audit report.
    (output / 'audit_summary.local.json').write_text(json.dumps({**summary, 'source_root': str(source.resolve())}, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({'files': len(records), 'code_or_docs': len(code), 'duplicate_groups': len(dupes)}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', required=True); p.add_argument('--output', required=True)
    a = p.parse_args(); audit(a.source, a.output)
