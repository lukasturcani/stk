import itertools as it
from functools import partial

import numpy as np
import pytest
import stk

from ..utilities import is_clone
from .check_atom_infos import check_atom_infos
from .check_atoms import check_atoms
from .check_bond_infos import check_bond_infos
from .check_bonds import check_bonds
from .check_building_block_counts import check_building_block_counts
from .check_building_blocks import check_building_blocks
from .check_edge_functional_groups import check_edge_functional_groups
from .check_edges import check_edges
from .check_position_matrix import check_position_matrix
from .check_vertex_edges import check_vertex_edges
from .check_vertices import check_vertices


@pytest.mark.skip
def test_with_placement_results(construction_state):
    building_blocks = tuple(
        map(
            construction_state.get_building_block,
            range(construction_state.get_num_vertices()),
        )
    )
    placement_results = tuple(
        map(
            partial(get_placement_result, construction_state),
            building_blocks,
        )
    )

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
    return stk.PlacementResult(
        position_matrix=np.array(
            [[i, i, i] for i in range(building_block.get_num_atoms())],
            dtype=np.float64,
        ),
        functional_group_edges={
            fg_id: next(edges)
            for fg_id in range(building_block.get_num_functional_groups())
        },
    )
