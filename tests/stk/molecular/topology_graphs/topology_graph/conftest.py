import numpy as np
import pytest
import stk


class VertexTestCase:
    def __init__(
        self,
        vertex,
        position,
        num_edges,
        edge_ids,
        cell,
        vertices,
        edges,
        building_block,
        building_block_centroid,
        alignments,
        functional_group_edges,
    ):
        self.vertex = vertex
        self.vertices = vertices
        self.edges = edges
        self.alignments = alignments
        self.functional_group_edges = functional_group_edges


@pytest.fixture(
    params=[

    ],
)
def make_vertex(request):
    return request.param


@pytest.fxutre(
    params=[

    ],
)
def make_edges(request):
    return request.param


@pytest.fixture(
    params=[0, 1, 0.5, 2],
)
def scale(request):
    return request.param


@pytest.fixture(
    params=[

    ],
)
def position(request):
    return np.array(request.param)


@pytest.fixture(
    params=[

    ],
)
def edge_ids(request):
    return list(request.param)


@pytest.fixture(
    params=[
        (0, 0, 0),
    ],
)
def cell(request):
    return tuple(request.param)
