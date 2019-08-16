import stk
import os
from os.path import join
import numpy as np


from ..._test_utilities import _test_dump_and_load


test_dir = 'host_guest_complex_topology_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def test_vertex():
    complex1 = stk.host_guest_complex.Complex()

    complex2 = stk.host_guest_complex.Complex(
        axis=[1, 0, 0],
        angle=np.pi/4,
        displacement=23
    )


def _create_host(amine2, amine2_alt1, aldehyde3):
    four_plus_six = stk.cage.FourPlusSix()
    return stk.ConstructedMolecule(
        building_blocks=[amine2, amine2_alt1, aldehyde3],
        topology_graph=four_plus_six,
        building_block_vertices={
            amine2: four_plus_six.vertices[5:],
            amine2_alt1: four_plus_six.vertices[4:5],
            aldehyde3: four_plus_six.vertices[:4]
        }
    )


def _test_construction(complex_, i, n):
    complex_.write(join(test_dir, f'complex_{i}.mol'))


def test_host_guest_complex(
    amine2,
    amine2_alt1,
    aldehyde3,
    chained_c60
):

    host = _create_host(amine2, amine2_alt1, aldehyde3)
    n = 3
    for i in range(n):
        complex_ = stk.ConstructedMolecule(
            building_blocks=[host, chained_c60],
            topology_graph=stk.host_guest_complex.Complex(
                axis=[1, 0, 0],
                angle=2*np.pi*i/n,
                displacement=[2**i, 0, 0]
            )
        )

        _test_construction(complex_, i, n)
        _test_dump_and_load(test_dir, complex_)
