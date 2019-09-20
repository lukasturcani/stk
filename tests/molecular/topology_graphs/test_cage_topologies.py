import stk
import os
from os.path import join
import numpy as np


from ..._test_utilities import _test_dump_and_load


test_dir = 'cage_topology_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def _alignment(vertex, building_block, edges):
    fg_position = building_block.get_centroid(
        atom_ids=building_block.func_groups[0].get_bonder_ids()
    )
    fg_vector = stk.normalize_vector(
        fg_position - vertex.get_position()
    )

    def inner(edge_id):
        edge_vector = (
            edges[edge_id].get_position() - vertex.get_position()
        )
        return fg_vector @ stk.normalize_vector(edge_vector)

    return inner


def _test_placement(vertex, bb, vertices, edges):
    vertex.place_building_block(bb, vertices, edges)
    assert np.allclose(
        a=bb.get_centroid(bb.get_bonder_ids()),
        b=vertex.get_position(),
        atol=1e-8,
    )
    aligned = max(
        vertex.get_edge_ids(),
        key=_alignment(vertex, bb, edges)
    )
    assert aligned is vertex._edge_ids[vertex.get_aligner_edge()]


def _test_assignment(vertex, bb, vertices, edges):
    assignments = vertex.assign_func_groups_to_edges(
        building_block=bb,
        vertices=vertices,
        edges=edges
    )
    assert (
        assignments[0] == vertex._edge_ids[vertex.get_aligner_edge()]
    )


def test_vertex(
    tmp_amine2,
    tmp_aldehyde3,
    tmp_aldehyde4,
    tmp_aldehyde5
):
    topology_graphs = (
        stk.cage.SixPlusEight(),
        stk.cage.OnePlusOne(),
        stk.cage.TwoPlusTwo(),
        stk.cage.FourPlusFour(),
        stk.cage.TwelvePlusThirty(),
        stk.cage.TwoPlusFour(),
        stk.cage.ThreePlusSix(),
        stk.cage.FourPlusEight(),
        stk.cage.FivePlusTen(),
        stk.cage.SixPlusTwelve(),
        stk.cage.EightPlusSixteen(),
        stk.cage.TenPlusTwenty(),
        stk.cage.TwoPlusThree(),
        stk.cage.FourPlusSix(),
        stk.cage.FourPlusSix2(),
        stk.cage.SixPlusNine(),
        stk.cage.EightPlusTwelve(),
        stk.cage.TwentyPlusThirty()
    )
    building_blocks = {
        2: tmp_amine2,
        3: tmp_aldehyde3,
        4: tmp_aldehyde4,
        5: tmp_aldehyde5
    }
    for topology_graph in topology_graphs:
        vertices = topology_graph.vertices
        edges = topology_graph.edges
        for vertex in topology_graph.vertices:
            bb = building_blocks[vertex.get_num_edges()]
            _test_placement(vertex, bb, vertices, edges)
            _test_assignment(vertex, bb, vertices, edges)


def test_topologies(
    tmp_six_plus_eight,
    tmp_one_plus_one,
    tmp_two_plus_two,
    tmp_four_plus_four,
    tmp_twelve_plus_thirty,
    tmp_two_plus_four,
    tmp_three_plus_six,
    tmp_four_plus_eight,
    tmp_five_plus_ten,
    tmp_six_plus_twelve,
    tmp_eight_plus_sixteen,
    tmp_ten_plus_twenty,
    tmp_two_plus_three,
    tmp_four_plus_six,
    tmp_four_plus_six2,
    tmp_six_plus_nine,
    tmp_eight_plus_twelve,
    tmp_twenty_plus_thirty
):
    cages = (
        (tmp_six_plus_eight, 6, 8),
        (tmp_one_plus_one, 1, 1),
        (tmp_two_plus_two, 2, 2),
        (tmp_four_plus_four, 4, 4),
        (tmp_twelve_plus_thirty, 12, 30),
        (tmp_two_plus_four, 2, 4),
        (tmp_three_plus_six, 3, 6),
        (tmp_four_plus_eight, 4, 8),
        (tmp_five_plus_ten, 5, 10),
        (tmp_six_plus_twelve, 6, 12),
        (tmp_eight_plus_sixteen, 8, 16),
        (tmp_ten_plus_twenty, 10, 20),
        (tmp_two_plus_three, 2, 3),
        (tmp_four_plus_six, 4, 6),
        (tmp_four_plus_six2, 4, 6),
        (tmp_six_plus_nine, 6, 9),
        (tmp_eight_plus_twelve, 8, 12),
        (tmp_twenty_plus_thirty, 20, 30),
    )
    for cage, num_expected_bb1s, num_expected_bb2s in cages:
        bb1, bb2 = sorted(
            cage.get_building_blocks(),
            key=lambda bb: len(bb.func_groups),
            reverse=True
        )
        num_expected_bbs = {
            bb1: num_expected_bb1s,
            bb2: num_expected_bb2s
        }
        _test_construction(cage, num_expected_bbs)
        _test_dump_and_load(test_dir, cage)


def test_alignments(amine2, amine2_alt3, aldehyde3, aldehyde3_alt3):
    building_blocks = [amine2, amine2_alt3, aldehyde3, aldehyde3_alt3]
    for fg in range(3):
        four_plus_six = stk.cage.FourPlusSix(
            vertex_alignments={3: fg},
        )
        c = stk.ConstructedMolecule(
            building_blocks=building_blocks,
            topology_graph=four_plus_six,
            building_block_vertices={
                amine2: four_plus_six.vertices[4:9],
                amine2_alt3: four_plus_six.vertices[9:],
                aldehyde3: four_plus_six.vertices[:3],
                aldehyde3_alt3: four_plus_six.vertices[3:4]
            }
        )
        c.write(join(test_dir, f'4p6_valignment_{fg}.mol'))

    four_plus_six = stk.cage.FourPlusSix(
        vertex_alignments={9: 1}
    )
    c = stk.ConstructedMolecule(
        building_blocks=building_blocks,
        topology_graph=four_plus_six,
        building_block_vertices={
            amine2: four_plus_six.vertices[4:9],
            amine2_alt3: four_plus_six.vertices[9:],
            aldehyde3: four_plus_six.vertices[:3],
            aldehyde3_alt3: four_plus_six.vertices[3:4]
        }
    )
    c.write(join(test_dir, f'4p6_edge_alignment.mol'))


def _test_construction(cage, num_expected_bbs):
    name = cage.topology_graph.__class__.__name__
    cage.write(join(test_dir, f'{name}_{len(num_expected_bbs)}.mol'))

    for bb in cage.get_building_blocks():
        assert cage.building_block_counter[bb] == num_expected_bbs[bb]
        # This test only holds true when each building block is
        # involved in every construction bond.
        if len(num_expected_bbs) < 3:
            assert (
                len(cage.construction_bonds) ==
                cage.building_block_counter[bb] * len(bb.func_groups)
            )
    num_deleters = sum(
        len(fg.deleters)*cage.building_block_counter[bb]
        for bb in cage.get_building_blocks() for fg in bb.func_groups
    )
    num_bb_atoms = sum(
        len(bb.atoms)*cage.building_block_counter[bb]
        for bb in cage.get_building_blocks()
    )
    num_bb_bonds = sum(
        len(bb.bonds)*cage.building_block_counter[bb]
        for bb in cage.get_building_blocks()
    )
    # Check that the correct number of bonds got made.
    assert (
        len(cage.construction_bonds) == len(cage.topology_graph.edges)
    )
    # Check correct total number of atoms.
    assert len(cage.atoms) == num_bb_atoms - num_deleters
    # Check correct total number of bonds.
    assert (
        len(cage.bonds) ==
        num_bb_bonds + len(cage.construction_bonds) - num_deleters
    )
    # Check window attributes got added.
    assert cage.num_windows == cage.topology_graph.num_windows
    assert (
        cage.num_window_types == cage.topology_graph.num_window_types
    )


def test_multi_bb(
    amine2,
    amine2_alt1,
    amine2_alt2,
    aldehyde3,
    aldehyde3_alt1,
    aldehyde3_alt2
):
    building_blocks = [
        amine2,
        amine2_alt1,
        amine2_alt2,
        aldehyde3,
        aldehyde3_alt1,
        aldehyde3_alt2
    ]

    four_plus_six = stk.cage.FourPlusSix()
    c = stk.ConstructedMolecule(
        building_blocks=building_blocks,
        topology_graph=four_plus_six,
        building_block_vertices={
            aldehyde3: four_plus_six.vertices[0:1],
            aldehyde3_alt1: four_plus_six.vertices[1:2],
            aldehyde3_alt2: four_plus_six.vertices[2:4],
            amine2: four_plus_six.vertices[4:6],
            amine2_alt1: four_plus_six.vertices[6:7],
            amine2_alt2: four_plus_six.vertices[7:]
        }
    )
    c.write(join(test_dir, 'multi_bb.mol'))
    num_expected_bbs = {
        amine2: 2,
        amine2_alt1: 1,
        amine2_alt2: 3,
        aldehyde3: 1,
        aldehyde3_alt1: 1,
        aldehyde3_alt2: 2
    }
    _test_construction(c, num_expected_bbs)
    _test_dump_and_load(test_dir, c)
