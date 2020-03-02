import pytest
import itertools as it
from functools import partial
import stk
import numpy as np

from ..utilities import is_clone
from .check_atoms import check_atoms
from .check_atom_infos import check_atom_infos
from .check_bonds import check_bonds
from .check_bond_infos import check_bond_infos
from .check_building_blocks import check_building_blocks
from .check_building_block_counts import check_building_block_counts
from .check_edges import check_edges
from .check_vertex_edges import check_vertex_edges
from .check_edge_functional_groups import check_edge_functional_groups
from .check_position_matrix import check_position_matrix
from .check_vertices import check_vertices


PlacementResult = (
    stk.molecular.topology_graphs.implementations
    .utilities._PlacementResult
)


def test_with_placement_results(construction_state):
    building_blocks = tuple(map(
        construction_state.get_building_block,
        range(construction_state.get_num_vertices()),
    ))
    placement_results = tuple(map(
        partial(get_placement_result, construction_state),
        building_blocks,
    ))

    # Clone for testing immutability.
    clone = construction_state.clone()
    new_state = construction_state.with_placement_results(
        building_blocks=building_blocks,
        results=placement_results,
    )
    old_state = construction_state
    check_atoms(old_state, new_state, building_blocks)
    check_atom_infos(old_state, new_state, building_blocks)
    check_bonds(old_state, new_state, building_blocks)
    check_bond_infos(old_state, new_state, building_blocks)
    check_building_blocks(old_state, new_state)
    check_building_block_counts(old_state, new_state, building_blocks)
    check_edges(old_state, new_state)
    check_vertex_edges(old_state, new_state)
    check_position_matrix(old_state, new_state, placement_results)
    check_vertices(old_state, new_state)
    check_edge_functional_groups(
        old_state=old_state,
        new_state=new_state,
        building_blocks=building_blocks,
        placement_results=placement_results,
    )
    is_clone(construction_state, clone)


def get_placement_result(construction_state, building_block):
    edges = it.cycle(range(construction_state.get_num_edges()))
    return PlacementResult(
        position_matrix=np.array(
            [[i, i, i] for i in range(building_block.get_num_atoms())],
            dtype=np.float64,
        ),
        functional_group_edges={
            fg_id: next(edges)
            for fg_id
            in range(building_block.get_num_functional_groups())
        },
    )


@pytest.fixture
def construction_state(building_block_vertices):
    return stk.ConstructionState(
        building_block_vertices=building_block_vertices,
        edges=get_edges(building_block_vertices),
        scale=1,
    )


def get_edges(building_block_vertices):
    vertices = (
        vertex
        for vertices in building_block_vertices.values()
        for vertex in vertices
    )
    vertex1 = next(vertices)
    for id_, vertex2 in enumerate(vertices):
        yield stk.Edge(id_, vertex1, vertex2)
        vertex1 = vertex2


@pytest.fixture
def building_block_vertices(building_block1, building_block2):
    return {
        building_block1: (
            stk.Vertex(0, [0, 0, 0]),
            stk.Vertex(1, [10, 0, 0]),
        ),
        building_block2: (
            stk.Vertex(2, [20, 0, 0]),
            stk.Vertex(3, [30, 0, 0]),
        ),
    }


@pytest.fixture(
    params=(
        stk.BuildingBlock('BrCC', [stk.BromoFactory()]),
    ),
)
def building_block1(request):
    return request.param


@pytest.fixture(
    params=(
        stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
    ),
)
def building_block2(request):
    return request.param
