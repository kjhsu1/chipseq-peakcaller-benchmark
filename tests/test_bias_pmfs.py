import sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.append(str(ROOT))
sys.argv = ['prog', '--output_fasta', 'dummy.fa']
from Kenta_Stuff.Scripts import copy_chip_seq as cs


def test_bias_pmfs():
    seq = 'ACGTACGTACGT'
    cs.k = 4
    length = len(seq) - cs.k + 1
    tf_bias = cs.build_tf_bias_pmf(length, [2, 6], sigma=1.0, enrichment=0.5)
    assert tf_bias.shape[0] == length
    np.testing.assert_almost_equal(tf_bias.sum(), 1.0)

    gc_bias = cs.build_gc_bias_pmf(seq, {'csv': 'data/gc_bias_example.csv'})
    assert gc_bias.shape[0] == length
    np.testing.assert_almost_equal(gc_bias.sum(), 1.0)

    acc_bias = cs.build_accessibility_bias_pmf(length, 'tests/data/example.bed', 2.0)
    assert acc_bias.shape[0] == length
    np.testing.assert_almost_equal(acc_bias.sum(), 1.0)

    base = np.ones(length)
    final = base * tf_bias * gc_bias * acc_bias
    final = final / final.sum()
    np.testing.assert_almost_equal(final.sum(), 1.0)
