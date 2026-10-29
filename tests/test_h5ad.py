import importlib.util
import tempfile
import unittest
import contextlib
import uuid
from pathlib import Path

import h5py
import numpy as np
from scipy import io, sparse

spec = importlib.util.spec_from_file_location('bridge', Path(__file__).resolve().parents[1] / 'scripts/h5ad_to_mtx.py')
bridge = importlib.util.module_from_spec(spec); spec.loader.exec_module(bridge)
TEST_TMP = Path(__file__).resolve().parents[1] / '.test-output' / 'python'
TEST_TMP.mkdir(parents=True, exist_ok=True)


@contextlib.contextmanager
def fixture_dir():
    # Default inherited ACLs work in Windows managed workspaces; mkdtemp's
    # private ACL may hide files from the subprocess sandbox token.
    path = TEST_TMP / str(uuid.uuid4())
    path.mkdir()
    yield str(path)


class H5adTest(unittest.TestCase):
    def fixture(self, path, encoding):
        x = np.array([[1, 0, 3], [0, 2, 4]])
        with h5py.File(path, 'w') as f:
            for name, ids in [('obs', ['c1', 'c2']), ('var', ['A', 'B', 'C'])]:
                g = f.create_group(name); g.attrs['_index'] = '_index'
                g.create_dataset('_index', data=np.array(ids, dtype='S'))
            cat = f['obs'].create_group('condition'); cat.create_dataset('codes', data=[0, 1]); cat.create_dataset('categories', data=np.array(['Control', 'Disease'], dtype='S'))
            layers = f.create_group('layers')
            if encoding == 'dense':
                layers.create_dataset('counts', data=x)
            else:
                m = sparse.csr_matrix(x) if encoding == 'csr_matrix' else sparse.csc_matrix(x)
                g = layers.create_group('counts'); g.attrs['shape'] = m.shape; g.attrs['encoding-type'] = encoding
                for key in ['data', 'indices', 'indptr']: g.create_dataset(key, data=getattr(m, key))
        return x

    def test_orientations(self):
        for enc in ['dense', 'csr_matrix', 'csc_matrix']:
            with self.subTest(enc=enc), fixture_dir() as tmp:
                p = Path(tmp); x = self.fixture(p / 'in.h5ad', enc)
                bridge.convert(p / 'in.h5ad', p / 'out')
                with (p / 'out/matrix.mtx').open('rb') as handle:
                    np.testing.assert_array_equal(io.mmread(handle).toarray(), x.T)
                self.assertIn('Disease', (p / 'out/metadata.tsv').read_text())

    def test_no_raw_count_fallback(self):
        with fixture_dir() as tmp:
            p = Path(tmp); self.fixture(p / 'in.h5ad', 'dense')
            with h5py.File(p / 'in.h5ad', 'a') as f: del f['layers/counts']
            with self.assertRaisesRegex(ValueError, 'Missing explicit raw counts'): bridge.convert(p / 'in.h5ad', p / 'out')

    def test_reject_noninteger(self):
        with fixture_dir() as tmp:
            p = Path(tmp); self.fixture(p / 'in.h5ad', 'dense')
            with h5py.File(p / 'in.h5ad', 'a') as f:
                del f['layers/counts']; f['layers'].create_dataset('counts', data=np.ones((2, 3)) * .5)
            with self.assertRaisesRegex(ValueError, 'noninteger'): bridge.convert(p / 'in.h5ad', p / 'out')


if __name__ == '__main__': unittest.main()
