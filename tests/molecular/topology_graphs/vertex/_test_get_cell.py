import numpy as np


def test_get_cell(case_data):
    _test_get_cell(case_data.vertex, case_data.cell)


def _test_get_cell(vertex, cell):
    assert np.all(np.equal(vertex.get_cell(), cell))
