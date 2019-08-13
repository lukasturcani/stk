import stk
import os
from collections import namedtuple
from os.path import join
import numpy as np

from ..._test_utilities import _test_dump_and_load


test_dir = 'cof_topology_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def test_place_linear_building_block(tmp_amine2):
    topology_graphs = (
        stk.cof.Honeycomb((2, 2, 1)),
        stk.cof.Hexagonal((2, 2, 1)),
        stk.cof.Square((2, 2, 1)),
        stk.cof.Kagome((2, 2, 1)),
        stk.cof.LinkerlessHoneycomb((2, 2, 1))
    )
    for topology_graph in topology_graphs:
        linear_vertex = next(
            v for v in topology_graph.vertices if len(v.edges) == 2
        )
        linear_vertex.place_building_block(tmp_amine2)
        bonder_centroid = tmp_amine2.get_centroid(
            atom_ids=tmp_amine2.get_bonder_ids()
        )
        assert np.allclose(
            a=bonder_centroid,
            b=linear_vertex.get_position(),
            atol=1e-6
        )


def test_place_nonlinear_building_block(
    tmp_aldehyde3,
    tmp_aldehyde4,
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
        3: tmp_aldehyde3,
        4: tmp_aldehyde4,
        6: tmp_aldehyde6
    }
    for topology_graph in topology_graphs:
        nonlinear_vertex = next(
            v for v in topology_graph.vertices if len(v.edges) > 2
        )
        bb = building_blocks[len(nonlinear_vertex.edges)]
        nonlinear_vertex.place_building_block(bb)
        bonder_centroid = tmp_aldehyde3.get_centroid(
            atom_ids=tmp_aldehyde3.get_bonder_ids()
        )
        assert np.allclose(
            a=bonder_centroid,
            b=nonlinear_vertex.get_position(),
            atol=1e-6
        )


def _test_construction(
    cof,
    num_expected_bbs,
    num_unreacted_fgs
):
    path = join(
        test_dir, f'{cof.topology_graph.__class__.__name__}.mol'
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


def test_alignments(amine2_alt3, aldehyde4_alt1):
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
    v0 = stk.cof.Kagome.vertices[0]
    vlast = stk.cof.Kagome.vertices[-1]
    for i in range(4):
        for periodic in (True, False):
            cof = stk.ConstructedMolecule(
                building_blocks=[amine2_alt3, aldehyde4_alt1],
                topology_graph=stk.cof.Kagome(
                    lattice_size=(3, 3, 1),
                    vertex_alignments={
                        v0: v0.edges[i],
                        vlast: vlast.edges[i % 2]
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
                num_unreacted_fgs=num_unreacted_fgs
            )
            _test_dump_and_load(test_dir, cof)


def test_multi_bb(
    amine2,
    amine2_alt1,
    amine2_alt2,
    amine2_alt3,
    aldehyde4,
    aldehyde4_alt1
):
    building_blocks = [
        amine2,
        amine2_alt1,
        amine2_alt2,
        amine2_alt3,
        aldehyde4,
        aldehyde4_alt1
    ]
    num_expected_bbs = {
        amine2: 1,
        amine2_alt1: 1,
        amine2_alt2: 2,
        amine2_alt3: 2,
        aldehyde4: 2,
        aldehyde4_alt1: 1
    }
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
        cof = stk.ConstructedMolecule(
            building_blocks=building_blocks,
            topology_graph=kagome,
            building_block_vertices={
                amine2: kagome.vertices[3:4],
                amine2_alt1: kagome.vertices[4:5],
                amine2_alt2: kagome.vertices[5:7],
                amine2_alt3: kagome.vertices[7:],
                aldehyde4: kagome.vertices[:2],
                aldehyde4_alt1: kagome.vertices[2:3]
            }
        )
        kind = '_periodic' if periodic else ''
        cof.write(join(test_dir, f'multi_bb{kind}.mol'))
        num_unreacted_fgs = (
            periodic_unreacted if periodic else island_unreacted
        )
        _test_construction(cof, num_expected_bbs, num_unreacted_fgs)
        _test_dump_and_load(test_dir, cof)


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
    tmp_periodic_linkerless_honeycomb
):

    COFData = namedtuple(
        'COFData',
        [
            'cof',
            'num_linkers',
            'num_building_blocks',
            'num_unreacted_linker_fgs',
            'num_unreacted_building_block_fgs'
            'periodic'
        ]
    )
    cofs = (
        COFData(tmp_honeycomb, 3*9, 2*9, 23, 23, False),
        COFData(tmp_periodic_honeycomb, 3*9, 2*9, 0, 0, True),
        COFData(tmp_kagome, 6*9, 3*9, 11, 11, False),
        COFData(tmp_periodic_kagome, 6*9, 3*9, 0, 0, True),
        COFData(tmp_hexagonal, 12*9, 4*9, 34, 34, False),
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
            num_unreacted_fgs=num_unreacted_fgs
        )
        _test_dump_and_load(cof.cof)
