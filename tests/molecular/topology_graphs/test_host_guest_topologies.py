import stk
import os
from os.path import join
import numpy as np


from ..._test_utilities import _test_dump_and_load, _compare_with_valid


test_dir = 'host_guest_topology_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def _test_placement(vertex, bb, vertices, edges):
    vertex.place_building_block(bb, vertices, edges)
    assert np.allclose(
        a=bb.get_centroid(),
        b=vertex.get_position(),
        atol=1e-5
    )


def _test_alignment(vertex, bb, target):
    atom1, atom2 = bb.get_atom_positions([0, 1])
    assert np.allclose(
        a=stk.normalize_vector(atom1 - atom2),
        b=target,
        atol=1e-8
    )


def test_vertex(tmp_four_plus_six, tmp_bromine2):
    host = tmp_four_plus_six
    guest = tmp_bromine2
    atom1, atom2 = guest.get_atom_positions([0, 1])
    guest_start = stk.normalize_vector(atom1 - atom2)
    complex2 = stk.host_guest.Complex(
        guest_start=guest_start,
        guest_target=[1, 1, 1],
        displacement=[23, 5, 21]
    )
    complexes = (
        (stk.host_guest.Complex(), guest_start),
        (complex2, stk.normalize_vector([1, 1, 1]))
    )
    for complex_, target in complexes:
        v1, v2 = complex_.vertices
        vertices = complex_.vertices
        edges = complex_.edges
        _test_placement(v1, host, vertices, edges)
        _test_placement(v2, guest, vertices, edges)
        _test_alignment(v2, guest, target)


def _create_host(amine2, amine2_alt3, aldehyde3):
    four_plus_six = stk.cage.FourPlusSix()
    return stk.ConstructedMolecule(
        building_blocks=[amine2, amine2_alt3, aldehyde3],
        topology_graph=four_plus_six,
        building_block_vertices={
            amine2: four_plus_six.vertices[5:],
            amine2_alt3: four_plus_six.vertices[4:5],
            aldehyde3: four_plus_six.vertices[:4]
        }
    )


def _test_construction(complex_, i, num_expected_bbs):
    complex_.write(join(test_dir, f'complex_{i}.mol'))

    assert len(complex_.building_block_counter) == 2
    for bb, num_expected in num_expected_bbs.items():
        assert complex_.building_block_counter[bb] == num_expected

    assert len(complex_.construction_bonds) == 0
    bb_atoms = sum(
        len(bb.atoms) for bb in complex_.get_building_blocks()
    )
    assert len(complex_.atoms) == bb_atoms
    bb_bonds = sum(
        len(bb.bonds) for bb in complex_.get_building_blocks()
    )
    assert len(complex_.bonds) == bb_bonds


def test_complex(
    amine2,
    amine2_alt3,
    aldehyde3,
    chained_c60,
    valid_host_guest_dir
):

    host = _create_host(amine2, amine2_alt3, aldehyde3)
    num_expected_bbs = {host: 1, chained_c60: 1}
    n = 3
    for i in range(n):
        x, y = np.cos(i*2*np.pi/n), np.sin(i*2*np.pi/n)
        complex_ = stk.ConstructedMolecule(
            building_blocks=[host, chained_c60],
            topology_graph=stk.host_guest.Complex(
                guest_start=chained_c60.get_direction(),
                guest_target=[x, y, 0],
                displacement=[3**i, 0, 0]
            )
        )
        _test_construction(complex_, i, num_expected_bbs)
        _test_dump_and_load(test_dir, complex_, str(i))
        _compare_with_valid(valid_host_guest_dir, complex_, str(i))
