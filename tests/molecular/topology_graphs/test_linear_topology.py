import os
import stk
import numpy as np
from scipy.spatial.distance import euclidean
from os.path import join

from ..._test_utilities import _test_dump_and_load


test_dir = 'linear_topology_tests_output'
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
            range(len(bb.func_groups)), key=_fg_distance(edge, bb)
        )
        assert assignments[closest] == edge.id


def test_vertex(tmp_amine2):
    chain = stk.polymer.Linear(
        repeating_unit='AB',
        num_repeating_units=3
    )

    for vertex in chain.vertices:
        _test_placement(vertex, tmp_amine2)
        _test_assignment(vertex, tmp_amine2)


def _test_construction(
    polymer,
    repeat_units,
    num_lost_bonds_per_join
):
    fg = (
        next(polymer.get_building_blocks()).func_groups[0].fg_type.name
    )
    polymer.write(join(test_dir, f'polymer_{fg}.mol'))

    for bb in polymer.get_building_blocks():
        assert (
            polymer.building_block_counter[bb] == repeat_units
        )

    monomer_joins = 2*repeat_units - 1
    bonds_per_join = len(bb.func_groups[0].bonders)
    assert (
        len(polymer.construction_bonds) == monomer_joins*bonds_per_join
    )

    deleters_per_join = sum(
        len(bb.func_groups[0].deleters)
        for bb in polymer.get_building_blocks()
    )
    num_bb_atoms = sum(
        len(bb.atoms) for bb in polymer.get_building_blocks()
    )
    expected_atoms = (
        num_bb_atoms*repeat_units - deleters_per_join*monomer_joins
    )
    assert len(polymer.atoms) == expected_atoms

    num_bb_bonds = sum(
        len(bb.bonds) for bb in polymer.get_building_blocks()
    )
    expected_bonds = (
        num_bb_bonds*repeat_units -
        num_lost_bonds_per_join*monomer_joins
    )
    assert len(polymer.bonds) == expected_bonds


def test_construction(amine2, aldehyde2, boronic_acid2, diol2):
    repeat_units = 3

    polymers = (
        stk.ConstructedMolecule(
            building_blocks=[amine2, aldehyde2],
            topology_graph=stk.polymer.Linear(
                repeating_unit='AB',
                orientations=[1, 1],
                num_repeating_units=repeat_units
            )
        ),
        stk.ConstructedMolecule(
            building_blocks=[boronic_acid2, diol2],
            topology_graph=stk.polymer.Linear(
                repeating_unit='AB',
                num_repeating_units=repeat_units
            )
        )

    )
    lost_bonds_per_join = (2, 4)
    polymer_data = zip(polymers, lost_bonds_per_join)
    for polymer, num_lost_bonds_per_join in polymer_data:
        _test_construction(
            polymer=polymer,
            repeat_units=repeat_units,
            num_lost_bonds_per_join=num_lost_bonds_per_join
        )
        _test_dump_and_load(test_dir, polymer)
