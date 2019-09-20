import os
import stk
import numpy as np
import rdkit.Chem.AllChem as rdkit
from os.path import join

from ..._test_utilities import _test_dump_and_load, _compare_with_valid


test_dir = 'rotaxane_topology_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def _test_axle_placement(vertex,  bb, vertices, edges):
    vertex.place_building_block(bb, vertices, edges)
    assert np.allclose(
        a=vertex.get_position(),
        b=bb.get_centroid(),
        atol=1e-8
    )


def _cycle_atoms(bb):
    rdkit_mol = bb.to_rdkit_mol()
    return max(rdkit.GetSymmSSSR(rdkit_mol), key=len)


def _test_cycle_placement(vertex, bb, vertices, edges):
    vertex.place_building_block(bb, vertices, edges)
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
    vertices = rotaxane.vertices
    edges = rotaxane.edges
    axle, *cycles = rotaxane.vertices
    _test_axle_placement(axle, tmp_polymer, vertices, edges)
    for vertex in cycles:
        _test_cycle_placement(vertex, tmp_macrocycle, vertices, edges)


def _test_construction(test_dir, filename, rotaxane_data):
    rotaxane = rotaxane_data.rotaxane
    num_expected_bbs = rotaxane_data.num_expected_bbs
    rotaxane.write(join(test_dir, filename))

    assert (
        len(rotaxane.building_block_counter) == len(num_expected_bbs)
    )
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


class _RotaxaneData:
    def __init__(
        self,
        rotaxane,
        num_expected_bbs
    ):
        self.rotaxane = rotaxane
        self.num_expected_bbs = num_expected_bbs


def test_construction(
    tmp_polymer,
    tmp_macrocycle,
    tmp_macrocycle_alt1,
    valid_rotaxane_dir

):
    rotaxanes = (
        _RotaxaneData(
            rotaxane=stk.ConstructedMolecule(
                building_blocks=[tmp_polymer, tmp_macrocycle],
                topology_graph=stk.rotaxane.NRotaxane('A', 5)
            ),
            num_expected_bbs={
                tmp_polymer: 1,
                tmp_macrocycle: 5
            }
        ),

        _RotaxaneData(
            rotaxane=stk.ConstructedMolecule(
                building_blocks=[
                    tmp_polymer,
                    tmp_macrocycle,
                    tmp_macrocycle_alt1
                ],
                topology_graph=stk.rotaxane.NRotaxane('AAB', 5)
            ),
            num_expected_bbs={
                tmp_polymer: 1,
                tmp_macrocycle: 10,
                tmp_macrocycle_alt1: 5
            }
        ),

        _RotaxaneData(
            rotaxane=stk.ConstructedMolecule(
                building_blocks=[tmp_polymer, tmp_macrocycle],
                topology_graph=stk.rotaxane.NRotaxane('A', 1)
            ),
            num_expected_bbs={
                tmp_polymer: 1,
                tmp_macrocycle: 1
            }
        ),
    )

    for i, rotaxane in enumerate(rotaxanes):
        i = str(i)
        _test_construction(test_dir, f'{i}.mol', rotaxane)
        _test_dump_and_load(test_dir, rotaxane.rotaxane, i)
        _compare_with_valid(valid_rotaxane_dir, rotaxane.rotaxane, i)
