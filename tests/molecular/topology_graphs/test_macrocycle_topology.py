import os
from os.path import join
import stk
import numpy as np
from scipy.spatial.distance import euclidean


from ..._test_utilities import _test_dump_and_load


test_dir = 'cyclic_topology_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def _test_placement(vertex, bb):
    vertex.place_building_block(bb)
    assert np.allclose(
        a=vertex.get_position(),
        b=bb.get_centroid(bb.get_bonder_ids()),
        atol=1e-8
    )


def _fg_distance(edge, bb):
    edge_position = edge.get_position()

    def inner(fg_id):
        fg = bb.func_groups[fg_id]
        fg_position = bb.get_centroid(fg.get_bonder_ids())
        return euclidean(edge_position, fg_position)

    return inner


def _test_assignment(vertex, bb):
    assignments = vertex.assign_func_groups_to_edges(bb)
    for edge in vertex.edges:
        closest = min(
            range(len(bb.func_groups)),
            key=_fg_distance(edge, bb)
        )
        assert assignments[closest] == edge.id


def test_vertex(tmp_amine2):
    cycle = stk.macrocycle.Macrocycle(
        repeating_unit='AB',
        num_repeating_units=3
    )
    for vertex in cycle.vertices:
        _test_placement(vertex, tmp_amine2)
        _test_assignment(vertex, tmp_amine2)


def _test_construction(tmp_macrocycle):
    repeat_units = 3
    tmp_macrocycle.write(join(test_dir, f'macrocycle.mol'))

    assert len(tmp_macrocycle.building_block_vertices) == 2
    for bb in tmp_macrocycle.get_building_blocks():
        assert (
            tmp_macrocycle.building_block_counter[bb] == repeat_units
        )

    monomer_joins = 2*repeat_units
    assert len(tmp_macrocycle.construction_bonds) == monomer_joins

    deleters_per_join = sum(
        len(bb.func_groups[0].deleters)
        for bb in tmp_macrocycle.get_building_blocks()
    )
    num_bb_atoms = sum(
        len(bb.atoms) for bb in tmp_macrocycle.get_building_blocks()
    )
    expected_atoms = (
        num_bb_atoms*repeat_units - deleters_per_join*monomer_joins
    )
    assert len(tmp_macrocycle.atoms) == expected_atoms

    num_bb_bonds = sum(
        len(bb.bonds) for bb in tmp_macrocycle.get_building_blocks()
    )
    expected_bonds = num_bb_bonds*repeat_units - monomer_joins
    assert len(tmp_macrocycle.bonds) == expected_bonds


def test_construction(tmp_macrocycle, tmp_macrocycle_alt1):
    macrocycles = (
        tmp_macrocycle,
        tmp_macrocycle_alt1
    )

    for macrocycle in macrocycles:
        _test_construction(macrocycle)
        _test_dump_and_load(test_dir, macrocycle)
