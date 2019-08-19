import os
import stk
import numpy as np
import rdkit.Chem.AllChem as rdkit
from os.path import join

from ..._test_utilities import _test_dump_and_load


test_dir = 'rotaxane_topology_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def _test_axle_placement(vertex,  bb):
    vertex.place_building_block(bb)
    assert np.allclose(
        a=vertex.get_position(),
        b=bb.get_centroid(),
        atol=1e-8
    )


def _cycle_atoms(bb):
    rdkit_mol = bb.to_rdkit_mol()
    return max(rdkit.GetSymmSSSR(rdkit_mol), key=len)


def _test_cycle_placement(vertex, bb):
    vertex.place_building_block(bb)
    cycle_atoms = _cycle_atoms(bb)
    assert np.allclose(
        a=vertex.get_position(),
        b=bb.get_centroid(cycle_atoms),
        atol=1e-8
    )
    assert np.allclose(
        a=[1, 0, 0],
        b=bb.get_plane_normal(cycle_atoms),
        atol=1e-8
    )


def test_vertex(tmp_polymer, tmp_macrocycle):
    rotaxane = stk.rotaxane.NRotaxane('A', 4)
    axle, *cycles = rotaxane.vertices
    _test_axle_placement(axle, tmp_polymer)
    for vertex in cycles:
        _test_cycle_placement(vertex, tmp_macrocycle)


def _test_construction(test_dir, num_expected_bbs, rotaxane):
    rotaxane.write(join(test_dir, 'rotaxane.mol'))

    assert len(rotaxane.building_block_counter) == 2
    for bb in num_expected_bbs:
        assert (
            rotaxane.building_block_counter[bb] == num_expected_bbs[bb]
        )

    assert len(rotaxane.construction_bonds) == 0
    num_bb_atoms = sum(
        len(bb.atoms)*rotaxane.building_block_counter[bb]
        for bb in rotaxane.get_building_blocks()
    )
    assert len(rotaxane.atoms) == num_bb_atoms
    num_bb_bonds = sum(
        len(bb.bonds)*rotaxane.building_block_counter[bb]
        for bb in rotaxane.get_building_blocks()
    )
    assert len(rotaxane.bonds) == num_bb_bonds


def test_construction(tmp_rotaxane):
    polymer = next(
        bb for bb in tmp_rotaxane.get_building_blocks()
        if isinstance(bb.topology_graph, stk.polymer.Linear)
    )
    cycle = next(
        bb for bb in tmp_rotaxane.get_building_blocks()
        if bb is not polymer
    )
    num_expected_bbs = {
        polymer: 1,
        cycle: 5
    }
    _test_construction(test_dir, num_expected_bbs, tmp_rotaxane)
    _test_dump_and_load(test_dir, tmp_rotaxane)
