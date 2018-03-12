from ..molecular import (StructUnit3, StructUnit2, Periodic,
                         Honeycomb, Hexagonal, Square, Kagome)
import os

test_dir = 'cof_topology_tests'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def test_honeycomb():
    bb1 = StructUnit2(
        '/home/lukas/databases/liverpool_refined/amines_2f/0.mol')
    bb2 = StructUnit3(
        '/home/lukas/databases/liverpool_refined/aldehydes_3f/0.mol')
    cof = Periodic([bb1, bb2], Honeycomb())
    cof.write(os.path.join(test_dir, 'honeycomb.sdf'))


def test_hexagonal():
    bb1 = StructUnit2(
        '/home/lukas/databases/liverpool_refined/amines_2f/0.mol')
    bb2 = StructUnit3('/home/lukas/Dropbox/workspace/6_aldehyde.mol')
    cof = Periodic([bb1, bb2], Hexagonal())
    cof.write(os.path.join(test_dir, 'hexagonal.sdf'))


def test_square():
    bb1 = StructUnit2(
        '/home/lukas/databases/liverpool_refined/amines_2f/0.mol')
    bb2 = StructUnit3(
        '/home/lukas/databases/liverpool_refined/aldehydes_4f/59.mol')
    cof = Periodic([bb1, bb2], Square())
    cof.write(os.path.join(test_dir, 'square.sdf'))


def test_kagome():
    bb1 = StructUnit2(
        '/home/lukas/databases/liverpool_refined/amines_2f/0.mol')
    bb2 = StructUnit3(
        '/home/lukas/databases/liverpool_refined/aldehydes_4f/59.mol')
    cof = Periodic([bb1, bb2], Kagome())
    cof.write(os.path.join(test_dir, 'kagome.sdf'))
