import pytest
import numpy as np
import stk

from ..._test_case import _TestCase

vertices = stk.molecular.topology_graphs.cof.honeycomb


@pytest.fixture
def cof3(init_at_shifted_center, id, aligner_edge, cell, shift_params):
    return _TestCase(
        vertex=init_at_shifted_center(
            id=id,
            vertices=shift_params.vertices,
            cell_shifts=shift_params.cell_shifts,
            lattice_constants=shift_params.lattice_constants,
            aligner_edge=aligner_edge,
            cell=cell,
        ),
        id=id,
        position=shift_params.position,
        cell=cell,
    )


@pytest.fixture(
    params=(
        vertices._LinearCofVertex.init_at_shifted_center,
        vertices._NonLinearCofVertex.init_at_shifted_center,
    ),
)
def init_at_shifted_center(request):
    return request.param


class ShiftParams:
    def __init__(
        self,
        vertices,
        cell_shifts,
        lattice_constants,
        position,
    ):
        self.vertices = vertices
        self.cell_shifts = cell_shifts
        self.lattice_constants = lattice_constants
        self.position = position


@pytest.fixture(
    params=(
        ShiftParams(
            vertices=(
                stk.Vertex(1, [-1, 1, 2]),
                stk.Vertex(2, [1, -1, -2]),
            ),
            cell_shifts=(
                np.array([0, 0, 0]),
                np.array([0, 0, 0]),
            ),
            lattice_constants=(
                np.array([1, 0, 0]),
                np.array([0, 1, 0]),
                np.array([0, 0, 1]),
            ),
            position=np.array([0., 0., 0.]),
        ),
        ShiftParams(
            vertices=(
                stk.Vertex(1, [-1, 1, 2]),
                stk.Vertex(2, [1, -1, 0.]),
            ),
            cell_shifts=(
                np.array([10, -10, 10]),
                np.array([-10, 10, -10]),
            ),
            lattice_constants=(
                np.array([1, 0, 0]),
                np.array([0, 1, 0]),
                np.array([0, 0, 1]),
            ),
            position=np.array([0., 0., 1.]),
        ),
    ),
)
def shift_params(request):
    return request.param
