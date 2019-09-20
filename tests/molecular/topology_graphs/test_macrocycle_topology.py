import os
from os.path import join
import stk
import numpy as np
from scipy.spatial.distance import euclidean


from ..._test_utilities import _test_dump_and_load


test_dir = 'cyclic_topology_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def _test_placement(vertex, bb, vertices, edges):
    vertex.place_building_block(bb, vertices, edges)
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


def _test_assignment(vertex, bb, vertices, edges):
    assignments = vertex.assign_func_groups_to_edges(
        building_block=bb,
        vertices=vertices,
        edges=edges
    )
    for edge_id in vertex.get_edge_ids():
        edge = edges[edge_id]
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
    vertices = cycle.vertices
    edges = cycle.edges
    for vertex in cycle.vertices:
        _test_placement(vertex, tmp_amine2, vertices, edges)
        _test_assignment(vertex, tmp_amine2, vertices, edges)


def _test_construction(macrocycle_data):
    macrocycle = macrocycle_data.macrocycle
    repeating_unit = macrocycle_data.repeating_unit
    num_repeating_units = macrocycle_data.num_repeating_units
    num_expected_bbs = macrocycle_data.num_expected_bbs
    num_atoms_lost_per_join = macrocycle_data.num_atoms_lost_per_join
    num_bonds_lost_per_join = macrocycle_data.num_bonds_lost_per_join

    macrocycle.write(join(test_dir, f'macrocycle.mol'))

    assert len(macrocycle.building_block_vertices) == 2
    for bb in macrocycle.get_building_blocks():
        expected = num_expected_bbs[bb]
        assert macrocycle.building_block_counter[bb] == expected

    monomer_joins = len(repeating_unit)*num_repeating_units
    assert len(macrocycle.construction_bonds) == monomer_joins

    num_bb_atoms = sum(
        len(bb.atoms)*num_expected_bbs[bb]
        for bb in macrocycle.get_building_blocks()
    )
    expected_atoms = (
        num_bb_atoms - num_atoms_lost_per_join*monomer_joins
    )
    assert len(macrocycle.atoms) == expected_atoms

    num_bb_bonds = sum(
        len(bb.bonds)*num_expected_bbs[bb]
        for bb in macrocycle.get_building_blocks()
    )
    expected_bonds = (
        num_bb_bonds - num_bonds_lost_per_join*monomer_joins
    )
    assert len(macrocycle.bonds) == expected_bonds


class _MacrocycleData:
    def __init__(
        self,
        macrocycle,
        repeating_unit,
        num_repeating_units,
        num_expected_bbs,
        num_atoms_lost_per_join,
        num_bonds_lost_per_join
    ):
        self.macrocycle = macrocycle
        self.repeating_unit = repeating_unit
        self.num_repeating_units = num_repeating_units
        self.num_expected_bbs = num_expected_bbs
        self.num_atoms_lost_per_join = num_atoms_lost_per_join
        self.num_bonds_lost_per_join = num_bonds_lost_per_join


def test_construction():
    bb1 = stk.BuildingBlock('BrCCBr', ['bromine'])
    bb2 = stk.BuildingBlock('BrCNCBr', ['bromine'])
    macrocycles = (
        _MacrocycleData(
            macrocycle=stk.ConstructedMolecule(
                building_blocks=[bb1, bb2],
                topology_graph=stk.macrocycle.Macrocycle('AB', 3)
            ),
            repeating_unit='AB',
            num_repeating_units=3,
            num_expected_bbs={bb1: 3, bb2: 3},
            num_atoms_lost_per_join=2,
            num_bonds_lost_per_join=1
        ),
        _MacrocycleData(
            macrocycle=stk.ConstructedMolecule(
                building_blocks=[bb1, bb2],
                topology_graph=stk.macrocycle.Macrocycle('AAB', 3),
            ),
            repeating_unit='AAB',
            num_repeating_units=3,
            num_expected_bbs={bb1: 6, bb2: 3},
            num_atoms_lost_per_join=2,
            num_bonds_lost_per_join=1
        )
    )

    for macrocycle in macrocycles:
        _test_construction(macrocycle)
        _test_dump_and_load(test_dir, macrocycle.macrocycle)
