"""Extract a deterministic small technical demo from GSE163577 course data."""
import argparse
import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import io, sparse


def main(source, output, config):
    source, output = Path(source).resolve(), Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    ann = pd.read_csv(source / '12.Seurat/04.cellAnn.txt', sep='\t').set_index('id')['labels']
    # These are distinct AD hippocampus donors, verified from local GEO series metadata.
    donors = {'Disease.GSM4982084': 'AD1', 'Disease.GSM4982089': 'AD2', 'Disease.GSM4982090': 'AD3', 'Disease.GSM4982093': 'AD4'}
    manifest, metadata = [], []
    for sample, donor in donors.items():
        folder = source / '12.Seurat' / sample
        bars = pd.read_csv(folder / 'barcodes.tsv.gz', header=None, sep='\t')[0].astype(str).values
        cells = np.asarray([sample + '_' + b for b in bars])
        labels = ann.reindex(cells)
        ids = np.where(labels.values == 'Endothelial cells')[0]
        rng = np.random.default_rng(417)
        ids = np.sort(rng.choice(ids, min(240, len(ids)), replace=False))
        with gzip.open(folder / 'matrix.mtx.gz', 'rb') as handle:
            matrix = io.mmread(handle).tocsr()[:, ids]
        dest = output / sample; dest.mkdir(exist_ok=True)
        with (dest / 'matrix.mtx').open('wb') as handle:
            io.mmwrite(handle, matrix)
        with gzip.open(folder / 'features.tsv.gz', 'rt', encoding='utf-8') as f:
            (dest / 'features.tsv').write_text(f.read(), encoding='utf-8')
        (dest / 'barcodes.tsv').write_text('\n'.join(bars[ids]) + '\n', encoding='utf-8')
        manifest.append(dict(sample_id=sample, path=sample, donor_id=donor, condition='Disease', cell_type='Endothelial cells'))
        metadata.extend(dict(cell_id=c, sample_id=sample, donor_id=donor, condition='Disease', cell_type='Endothelial cells') for c in cells[ids])
    pd.DataFrame(manifest).to_csv(output / 'samples.tsv', sep='\t', index=False)
    pd.DataFrame(metadata).to_csv(output / 'metadata.tsv', sep='\t', index=False)
    base = Path(__file__).resolve().parents[1]
    cfg = yaml.safe_load((base / 'configs/fabp1_ppara_human.yaml').read_text(encoding='utf-8'))
    cfg['input'] = dict(type='10x', samples_file=str(output / 'samples.tsv'), metadata=str(output / 'metadata.tsv'), gene_column=2)
    cfg['analysis'].update(target_gene='GSN', readout_genes=['DLC1', 'PPARA'], cell_types=['Endothelial cells'], conditions=['Disease'], oe_factors=[2, 5, 10], technical_demo=True)
    cfg['network'].update(n_hvg=100, n_net=3, n_cells=120, td_max_iter=100, max_cells_per_donor=160, n_control_genes=1)
    cfg['qc'].update(min_genes=100, min_counts=500, max_mt_fraction=0.10, min_cells_network=100)
    cfg['enrichment'].update(go_gmt=str(source / '29.Enrich/../21.GSEA/c5.go.Hs.symbols.gmt'), kegg_gmt=str(source / '21.GSEA/c2.cp.kegg.Hs.symbols.gmt'), min_set_size=3)
    cfg['output']['directory'] = str(output.parent / 'course_results')
    Path(config).parent.mkdir(parents=True, exist_ok=True)
    Path(config).write_text(yaml.safe_dump(cfg, allow_unicode=True, sort_keys=False), encoding='utf-8')
    (output / 'provenance.json').write_text(json.dumps({'dataset': 'GSE163577', 'title': 'The human brain vasculature at single-cell resolution', 'technical_demo': True, 'selected_cell_type': 'Endothelial cells', 'cells': len(metadata), 'donors': donors, 'selection_seed': 417, 'labels': 'course-provided labels, not independent liver annotation', 'biological_scope': 'No liver fibrosis inference'}, ensure_ascii=False, indent=2), encoding='utf-8')
    print(f'Prepared {len(metadata)} cells from {len(donors)} AD donors; technical demo only.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', required=True); p.add_argument('--output', required=True); p.add_argument('--config', required=True)
    a = p.parse_args(); main(a.source, a.output, a.config)
