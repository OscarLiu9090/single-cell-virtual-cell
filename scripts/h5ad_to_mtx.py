"""Read only an explicitly named raw-count layer; export genes x cells.

Dependencies: h5py, scipy, pandas, numpy. No unpickling or AnnData installation.
"""
import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import io, sparse


def strings(values):
    return np.asarray([v.decode('utf-8') if isinstance(v, bytes) else str(v) for v in values])


def read_column(node):
    if isinstance(node, h5py.Dataset):
        a = node[()]
        return strings(a) if a.dtype.kind in 'SO' else a
    if 'codes' in node and 'categories' in node:
        codes = node['codes'][()]
        cats = read_column(node['categories'])
        return np.asarray([cats[c] if c >= 0 else None for c in codes], dtype=object)
    if 'values' in node and 'mask' in node:
        a = read_column(node['values']).astype(object)
        a[node['mask'][()].astype(bool)] = None
        return a
    raise ValueError(f'Unsupported AnnData column encoding: {node.name}')


def read_frame(group):
    key = group.attrs.get('_index', '_index')
    if isinstance(key, bytes):
        key = key.decode()
    index = strings(read_column(group[key]))
    columns = {name: read_column(group[name]) for name in group if name != key}
    return pd.DataFrame(columns, index=index)


def convert(path, output, layer='counts', gene_symbol_column=None):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise ValueError('Output directory must be empty; refusing to overwrite.')
    with h5py.File(path, 'r') as f:
        obs, var = read_frame(f['obs']), read_frame(f['var'])
        key = f'layers/{layer}'
        if key not in f:
            raise ValueError(f'Missing explicit raw counts layer {key}; X/raw.X is never assumed to be raw counts.')
        node = f[key]
        if isinstance(node, h5py.Dataset):
            # Dense datasets are imported row blocks to bound peak memory.
            x = sparse.vstack([sparse.csr_matrix(node[i:i + 1024]) for i in range(0, node.shape[0], 1024)], format='csr')
        else:
            encoding = node.attrs.get('encoding-type', '')
            if isinstance(encoding, bytes):
                encoding = encoding.decode()
            cls = {'csr_matrix': sparse.csr_matrix, 'csc_matrix': sparse.csc_matrix}.get(encoding)
            if cls is None:
                raise ValueError(f'Unsupported matrix encoding {encoding}')
            x = cls((node['data'][()], node['indices'][()], node['indptr'][()]), shape=tuple(node.attrs['shape']))
        if x.shape != (len(obs), len(var)):
            raise ValueError('AnnData matrix shape does not match cells x genes metadata.')
        if obs.index.has_duplicates or obs.index.isna().any():
            raise ValueError('Cell IDs must be unique and nonmissing.')
        if not np.isfinite(x.data).all() or (x.data < 0).any() or not np.allclose(x.data, np.rint(x.data), atol=1e-6, rtol=0):
            raise ValueError('Explicit counts layer contains invalid/noninteger values.')
        genes = var[gene_symbol_column].astype(str).values if gene_symbol_column else var.index.values
        if any(str(g) in ('', 'nan', 'None') for g in genes):
            raise ValueError('Missing gene symbols.')
        if any('\n' in str(v) or '\t' in str(v) for v in list(genes) + list(obs.index)):
            raise ValueError('Identifiers cannot contain tabs/newlines.')
        with (output / 'matrix.mtx').open('wb') as handle:
            io.mmwrite(handle, x.T.tocoo())
        (output / 'genes.txt').write_text('\n'.join(map(str, genes)) + '\n', encoding='utf-8')
        (output / 'cells.txt').write_text('\n'.join(obs.index) + '\n', encoding='utf-8')
        obs['cell_id'] = obs.index
        obs.to_csv(output / 'metadata.tsv', sep='\t', index=False)
        (output / 'conversion.json').write_text(json.dumps({'input': str(Path(path).resolve()), 'layer': key, 'input_shape': x.shape, 'output_orientation': 'genes_x_cells', 'gene_symbol_column': gene_symbol_column}, indent=2), encoding='utf-8')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--layer', default='counts')
    parser.add_argument('--gene-symbol-column')
    a = parser.parse_args()
    convert(a.input, a.output, a.layer, a.gene_symbol_column)
