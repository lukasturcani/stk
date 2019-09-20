import stk
import os
from collections import namedtuple
from os.path import join
import numpy as np

from ..._test_utilities import _test_dump_and_load, _compare_with_valid


test_dir = 'cof_topology_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def _alignment(vertex, building_block, vertices, edges):
    fg_position = building_block.get_centroid(
        atom_ids=building_block.func_groups[0].get_bonder_ids()
    )
    v1 = stk.normalize_vector(fg_position - vertex.get_position())

    def inner(edge_id):
        edge_position = edges[edge_id].get_position(vertex, vertices)
        v2 = edge_position - vertex.get_position()
        return v1 @ stk.normalize_vector(v2)

    return inner


def _test_placement(vertex, bb, vertices, edges):
    vertex.place_building_block(bb, vertices, edges)
    assert np.allclose(
        a=bb.get_centroid(bb.get_bonder_ids()),
        b=vertex.get_position(),
        atol=1e-8
    )
    aligned = max(
        vertex.get_edge_ids(),
        key=_alignment(vertex, bb, vertices, edges)
    )
    vertex_edges = list(vertex.get_edge_ids())
    assert aligned == vertex_edges[vertex.get_aligner_edge()]


def _angle(bb, edge, vertex, vertices):
    edge_vector = (
        edge.get_position(vertex, vertices) -
        bb.get_centroid(bb.get_bonder_ids())
    )

    def inner(fg_id):
        fg = bb.func_groups[fg_id]
        fg_vector = (
            bb.get_centroid(fg.get_bonder_ids()) -
            bb.get_centroid(bb.get_bonder_ids())
        )
        return stk.vector_angle(fg_vector, edge_vector)

    return inner


def _test_assignment(vertex, bb, vertices, edges):
    assignments = (
        vertex.assign_func_groups_to_edges(bb, vertices, edges)
    )
    vertex_edges = list(vertex.get_edge_ids())
    assert assignments[0] == vertex_edges[vertex.get_aligner_edge()]
    for edge_id in vertex.get_edge_ids():
        closest = min(
            range(len(bb.func_groups)),
            key=_angle(bb, edges[edge_id], vertex, vertices)
        )
        assert assignments[closest] == edge_id

    if len(bb.func_groups) == 2:
        not_aligner = next(
            e for e in vertex.get_edge_ids()
            if e != vertex_edges[vertex.get_aligner_edge()]
        )
        assert assignments[1] == not_aligner


def test_vertex(
    tmp_amine2,
    tmp_aldehyde3,
    tmp_aldehyde4_alt2,
    tmp_aldehyde6
):
    topology_graphs = (
        stk.cof.Honeycomb((2, 2, 1)),
        stk.cof.Hexagonal((2, 2, 1)),
        stk.cof.Square((2, 2, 1)),
        stk.cof.Kagome((2, 2, 1)),
        stk.cof.LinkerlessHoneycomb((2, 2, 1))
    )
    building_blocks = {
        2: tmp_amine2,
        3: tmp_aldehyde3,
        4: tmp_aldehyde4_alt2,
        6: tmp_aldehyde6
    }
    for topology_graph in topology_graphs:
        vertices = topology_graph.vertices
        edges = topology_graph.edges

        for vertex in topology_graph.vertices:
            bb = building_blocks[vertex.get_num_edges()]
            _test_placement(vertex, bb, vertices, edges)
            _test_assignment(vertex, bb, vertices, edges)


def _test_construction(
    cof,
    num_expected_bbs,
    num_unreacted_fgs,
    periodic
):
    is_periodic = '_periodic' if periodic else ''
    path = join(
        test_dir,
        f'{cof.topology_graph.__class__.__name__}{is_periodic}.mol'
    )
    cof.write(path)

    for bb in cof.get_building_blocks():
        assert cof.building_block_counter[bb] == num_expected_bbs[bb]
        # This test only holds true when each building block is
        # involved in every construction bond and the cof is
        # periodic.
        is_periodic = all(
            num == 0 for num in num_unreacted_fgs.values()
        )
        if len(num_expected_bbs) < 3 and is_periodic:
            assert (
                len(cof.construction_bonds) ==
                cof.building_block_counter[bb] * len(bb.func_groups)
            )

    # Check that the correct number of bonds got made.
    assert (
        len(cof.construction_bonds) ==
        # For every 2 unreacted functional groups there should be one
        # less construction bond.
        len(cof.topology_graph.edges)-sum(num_unreacted_fgs.values())/2
    )
    # Check correct total number of atoms.
    num_deleters = sum(
        len(fg.deleters)*cof.building_block_counter[bb]
        for bb in cof.get_building_blocks()
        for fg in bb.func_groups
    )
    # Remove the deleters from all the uncreated functional groups.
    # This assumes that all functional groups in a building block are
    # the same.
    num_deleters -= sum(
        len(bb.func_groups[0].deleters)
        for bb in num_unreacted_fgs
        for i in range(num_unreacted_fgs[bb])
    )
    num_bb_atoms = sum(
        len(bb.atoms)*cof.building_block_counter[bb]
        for bb in cof.get_building_blocks()
    )
    assert len(cof.atoms) == num_bb_atoms - num_deleters
    # Check correct total number of bonds.
    num_bb_bonds = sum(
        len(bb.bonds)*cof.building_block_counter[bb]
        for bb in cof.get_building_blocks()
    )
    assert (
        len(cof.bonds) ==
        num_bb_bonds + len(cof.construction_bonds) - num_deleters
    )


def test_alignments(amine2_alt3, aldehyde4_alt1, valid_cof_dir):
    num_expected_bbs = {
        amine2_alt3: 6*9,
        aldehyde4_alt1: 3*9
    }
    periodic_unreacted = {
        amine2_alt3: 0,
        aldehyde4_alt1: 0
    }
    island_unreacted = {
        amine2_alt3: 11,
        aldehyde4_alt1: 11
    }
    for i in range(4):
        for periodic in (True, False):
            cof = stk.ConstructedMolecule(
                building_blocks=[amine2_alt3, aldehyde4_alt1],
                topology_graph=stk.cof.Kagome(
                    lattice_size=(3, 3, 1),
                    vertex_alignments={
                        0: i,
                        len(stk.cof.Kagome.vertex_data)-1: i % 2
                    },
                    periodic=periodic
                )
            )
            kind = '_periodic' if periodic else ''
            cof.write(
                join(test_dir, f'aligning_{i}_{i%2}{kind}.mol')
            )
            num_unreacted_fgs = (
                periodic_unreacted if periodic else island_unreacted
            )
            _test_construction(
                cof=cof,
                num_expected_bbs=num_expected_bbs,
                num_unreacted_fgs=num_unreacted_fgs,
                periodic=periodic
            )
            _test_dump_and_load(
                test_dir=test_dir,
                mol=cof,
                name=f'aligning_{i}_{i%2}{kind}'
            )
            _compare_with_valid(
                test_dir=valid_cof_dir,
                mol=cof,
                name=f'aligning_{i}_{i%2}{kind}'
            )


def test_multi_bb(
    amine2,
    amine2_alt1,
    amine2_alt2,
    amine2_alt3,
    aldehyde4,
    aldehyde4_alt1,
    valid_cof_dir
):
    building_blocks = [
        amine2,
        amine2_alt1,
        amine2_alt2,
        amine2_alt3,
        aldehyde4,
        aldehyde4_alt1
    ]
    periodic_unreacted = {bb: 0 for bb in building_blocks}
    island_unreacted = {
        amine2: 0,
        amine2_alt1: 0,
        amine2_alt2: 3,
        amine2_alt3: 8,
        aldehyde4: 3,
        aldehyde4_alt1: 8
    }
    for periodic in (True, False):
        kagome = stk.cof.Kagome((3, 3, 1), periodic)
        di_verts = [
            v for v in kagome.vertices if v.get_num_edges() == 2
        ]
        tetra_verts = [
            v for v in kagome.vertices if v.get_num_edges() == 4
        ]
        cof = stk.ConstructedMolecule(
            building_blocks=building_blocks,
            topology_graph=kagome,
            building_block_vertices={
                amine2: di_verts[:4],
                amine2_alt1: di_verts[4:5],
                amine2_alt2: di_verts[5:7],
                amine2_alt3: di_verts[7:],
                aldehyde4: tetra_verts[:2],
                aldehyde4_alt1: tetra_verts[2:]
            }
        )
        num_expected_bbs = {
            amine2: len(di_verts[:4]),
            amine2_alt1: len(di_verts[4:5]),
            amine2_alt2: len(di_verts[5:7]),
            amine2_alt3: len(di_verts[7:]),
            aldehyde4: len(tetra_verts[:2]),
            aldehyde4_alt1: len(tetra_verts[2:])
        }
        kind = '_periodic' if periodic else ''
        cof.write(join(test_dir, f'multi_bb{kind}.mol'))
        num_unreacted_fgs = (
            periodic_unreacted if periodic else island_unreacted
        )
        _test_construction(
            cof=cof,
            num_expected_bbs=num_expected_bbs,
            num_unreacted_fgs=num_unreacted_fgs,
            periodic=periodic
        )
        _test_dump_and_load(test_dir, cof, 'multi')
        _compare_with_valid(valid_cof_dir, cof, 'multi')


def test_topologies(
    tmp_honeycomb,
    tmp_periodic_honeycomb,
    tmp_kagome,
    tmp_periodic_kagome,
    tmp_hexagonal,
    tmp_periodic_hexagonal,
    tmp_square,
    tmp_periodic_square,
    tmp_linkerless_honeycomb,
    tmp_periodic_linkerless_honeycomb,
    valid_cof_dir
):

    COFData = namedtuple(
        'COFData',
        [
            'cof',
            'num_linkers',
            'num_building_blocks',
            'num_unreacted_linker_fgs',
            'num_unreacted_building_block_fgs',
            'periodic'
        ]
    )
    cofs = (
        COFData(tmp_honeycomb, 3*9, 2*9, 6, 6, False),
        COFData(tmp_periodic_honeycomb, 3*9, 2*9, 0, 0, True),
        COFData(tmp_kagome, 6*9, 3*9, 11, 11, False),
        COFData(tmp_periodic_kagome, 6*9, 3*9, 0, 0, True),
        COFData(tmp_hexagonal, 12*9, 4*9, 23, 23, False),
        COFData(tmp_periodic_hexagonal, 12*9, 4*9, 0, 0, True),
        COFData(tmp_square, 2*9, 1*9, 6, 6, False),
        COFData(tmp_periodic_square, 2*9, 1*9, 0, 0, True),
        COFData(tmp_linkerless_honeycomb, 1*9, 1*9, 6, 6, False),
        COFData(tmp_periodic_linkerless_honeycomb, 9, 9, 0, 0, True)
    )

    for cof in cofs:
        linker, building_block = sorted(
            cof.cof.get_building_blocks(),
            key=lambda bb: len(bb.func_groups)
        )
        num_expected_bbs = {
            linker: cof.num_linkers,
            building_block: cof.num_building_blocks
        }
        num_unreacted_fgs = {
            linker: cof.num_unreacted_linker_fgs,
            building_block: cof.num_unreacted_building_block_fgs
        }
        _test_construction(
            cof=cof.cof,
            num_expected_bbs=num_expected_bbs,
            num_unreacted_fgs=num_unreacted_fgs,
            periodic=cof.periodic
        )
        _test_dump_and_load(test_dir, cof.cof)
        _compare_with_valid(valid_cof_dir, cof.cof)
