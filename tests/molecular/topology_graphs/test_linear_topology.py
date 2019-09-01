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


def _test_construction(polymer_data):
    polymer = polymer_data.polymer
    num_expected_bbs = polymer_data.num_expected_bbs
    num_repeating_units = polymer_data.num_repeating_units
    num_atoms_lost_per_join = polymer_data.num_atoms_lost_per_join
    num_bonds_lost_per_join = polymer_data.num_bonds_lost_per_join

    fg = (
        next(polymer.get_building_blocks()).func_groups[0].fg_type.name
    )
    polymer.write(join(test_dir, f'polymer_{fg}.mol'))

    for bb in polymer.get_building_blocks():
        assert (
            polymer.building_block_counter[bb] == num_expected_bbs[bb]
        )

    repeat_unit_size = len(polymer.topology_graph._repeating_unit)
    monomer_joins = repeat_unit_size*num_repeating_units - 1
    bonds_per_join = len(bb.func_groups[0].bonders)
    assert (
        len(polymer.construction_bonds) == monomer_joins*bonds_per_join
    )

    num_bb_atoms = sum(
        len(bb.atoms)*num_expected_bbs[bb]
        for bb in polymer.get_building_blocks()
    )
    expected_atoms = (
        num_bb_atoms - num_atoms_lost_per_join*monomer_joins
    )
    assert len(polymer.atoms) == expected_atoms

    num_bb_bonds = sum(
        len(bb.bonds)*num_expected_bbs[bb]
        for bb in polymer.get_building_blocks()
    )
    expected_bonds = (
        num_bb_bonds - num_bonds_lost_per_join*monomer_joins
    )
    assert len(polymer.bonds) == expected_bonds


class _PolymerData:
    def __init__(
        self,
        polymer,
        num_repeating_units,
        num_expected_bbs,
        num_atoms_lost_per_join,
        num_bonds_lost_per_join
    ):
        self.polymer = polymer
        self.num_repeating_units = num_repeating_units
        self.num_expected_bbs = num_expected_bbs
        self.num_atoms_lost_per_join = num_atoms_lost_per_join
        self.num_bonds_lost_per_join = num_bonds_lost_per_join


def test_construction(
    amine2,
    aldehyde2,
    boronic_acid2,
    diol2,
    tmp_bromine2,
    tmp_bromine2_alt1
):
    num_repeating_units = 3

    bb1 = stk.BuildingBlock('BrCCC', ['bromine'])
    bb2 = stk.BuildingBlock('BrCCCBr', ['bromine'])
    bb3 = stk.BuildingBlock('BrCCN', ['bromine'])

    polymers = (
        _PolymerData(
            polymer=stk.ConstructedMolecule(
                building_blocks=[amine2, aldehyde2],
                topology_graph=stk.polymer.Linear(
                    repeating_unit='AB',
                    orientations=[1, 1],
                    num_repeating_units=num_repeating_units
                )
            ),
            num_repeating_units=num_repeating_units,
            num_expected_bbs={
                amine2: num_repeating_units,
                aldehyde2: num_repeating_units
            },
            num_atoms_lost_per_join=3,
            num_bonds_lost_per_join=2
        ),

        _PolymerData(
            polymer=stk.ConstructedMolecule(
                building_blocks=[boronic_acid2, diol2],
                topology_graph=stk.polymer.Linear(
                    repeating_unit='AB',
                    num_repeating_units=num_repeating_units
                )
            ),
            num_repeating_units=num_repeating_units,
            num_expected_bbs={
                boronic_acid2: num_repeating_units,
                diol2: num_repeating_units
            },
            num_atoms_lost_per_join=6,
            num_bonds_lost_per_join=4
        ),

        _PolymerData(
            polymer=stk.ConstructedMolecule(
                building_blocks=[
                    tmp_bromine2,
                    tmp_bromine2_alt1
                ],
                topology_graph=stk.polymer.Linear(
                    repeating_unit='AAB',
                    num_repeating_units=num_repeating_units
                )
            ),
            num_repeating_units=num_repeating_units,
            num_expected_bbs={
                tmp_bromine2: 2*num_repeating_units,
                tmp_bromine2_alt1: num_repeating_units
            },
            num_atoms_lost_per_join=2,
            num_bonds_lost_per_join=1
        ),

        _PolymerData(
            polymer=stk.ConstructedMolecule(
                building_blocks=[bb1, bb2, bb3],
                topology_graph=stk.polymer.Linear('ABC', 1)
            ),
            num_repeating_units=1,
            num_expected_bbs={bb1: 1, bb2: 1, bb3: 1},
            num_atoms_lost_per_join=2,
            num_bonds_lost_per_join=1
        ),

    )

    for polymer_data in polymers:
        _test_construction(polymer_data)
        _test_dump_and_load(test_dir, polymer_data.polymer)
