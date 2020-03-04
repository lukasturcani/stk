import numpy as np


def test_get_cell(test_case):
    _test_get_cell(test_case.vertex, test_case.cell)


def _test_get_cell(vertex, cell):
    assert np.all(np.equal(vertex.get_cell(), cell))
