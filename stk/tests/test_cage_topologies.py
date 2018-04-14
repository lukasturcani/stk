from ..molecular.topologies.cage import *
from ..molecular import StructUnit3, StructUnit2, Cage
import os
from os.path import join

test_dir = 'cage_topology_tests'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)
data_dir = os.path.join(os.getcwd(), 'data', 'cage_topologies')


# 3 + 4 topology tests.
def test_SixPlusEight():
    bb1 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    bb2 = StructUnit3(join(data_dir, 'amine4.mol'))
    c = Cage([bb1, bb2], SixPlusEight())
    c.write(join(test_dir, 'SixPlusEight.mol'))


# 2 + 4 topolgy tests.
def test_TwoPlusFour():
    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit3(join(data_dir, 'amine4.mol'))
    c = Cage([bb1, bb2], TwoPlusFour())
    c.write(join(test_dir, 'TwoPlusFour.mol'))


def test_ThreePlusSix():
    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit3(join(data_dir, 'amine4.mol'))
    c = Cage([bb1, bb2], ThreePlusSix())
    c.write(join(test_dir, 'ThreePlusSix.mol'))


def test_FourPlusEight():
    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit3(join(data_dir, 'amine4.mol'))
    c = Cage([bb1, bb2], FourPlusEight())
    c.write(join(test_dir, 'FourPlusEight.mol'))


def test_FivePlusTen():
    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit3(join(data_dir, 'amine4.mol'))
    c = Cage([bb1, bb2], FivePlusTen())
    c.write(join(test_dir, 'FivePlusTen.mol'))


def test_SixPlusTwelve():
    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit3(join(data_dir, 'amine4.mol'))
    c = Cage([bb1, bb2], SixPlusTwelve())
    c.write(join(test_dir, 'SixPlusTwelve.mol'))


def test_EightPlusSixteen():
    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit3(join(data_dir, 'amine4.mol'))
    c = Cage([bb1, bb2], EightPlusSixteen())
    c.write(join(test_dir, 'EightPlusSixteen.mol'))


def test_TenPlusTwenty():
    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit3(join(data_dir, 'amine4.mol'))
    c = Cage([bb1, bb2], TenPlusTwenty())
    c.write(join(test_dir, 'TenPlusTwenty.mol'))


# 3 + 3 cage topologies.
def test_OnePlusOne():
    bb1 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    bb2 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = Cage([bb1, bb2], OnePlusOne())
    c.write(join(test_dir, 'OnePlusOne.mol'))


def test_TwoPlusTwo():
    bb1 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    bb2 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = Cage([bb1, bb2], TwoPlusTwo())
    c.write(join(test_dir, 'TwoPlusTwo.mol'))


def test_FourPlusFour():
    bb1 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    bb2 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = Cage([bb1, bb2], FourPlusFour())
    c.write(join(test_dir, 'FourPlusFour.mol'))


# 2 + 3 cage topologies.
def test_TwoPlusThree():
    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = Cage([bb1, bb2], TwoPlusThree())
    c.write(join(test_dir, 'TwoPlusThree.mol'))


def test_FourPlusSix():
    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = Cage([bb1, bb2], FourPlusSix())
    c.write(join(test_dir, 'FourPlusSix.mol'))


def test_multiFourPlusSix():
    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit2(join(data_dir, 'amine2_1.mol'))
    bb3 = StructUnit2(join(data_dir, 'amine2_2.mol'))

    bb4 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    bb5 = StructUnit3(join(data_dir, 'aldehyde3_1.mol'))
    bb6 = StructUnit3(join(data_dir, 'aldehyde3_2.mol'))

    c1 = Cage([bb1, bb2, bb3, bb4, bb5, bb6], FourPlusSix())

    c2 = Cage([bb1, bb2, bb3, bb4, bb5, bb6],
              FourPlusSix(bb_assignments={0: [0, 1],
                                          1: [2, 3, 4],
                                          2: [5],
                                          3: [0],
                                          4: [1, 2],
                                          5: [3]}))

    c1.write(join(test_dir, 'multi_FourPlusSix_1.mol'))
    c2.write(join(test_dir, 'multi_FourPlusSix_2.mol'))


def test_FourPlusSix2():
    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = Cage([bb1, bb2], FourPlusSix2())
    c.write(join(test_dir, 'FourPlusSix2.mol'))


def test_SixPlusNine():
    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = Cage([bb1, bb2], SixPlusNine())
    c.write(join(test_dir, 'SixPlusNine.mol'))


def test_EightPlusTwelve():
    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = Cage([bb1, bb2], EightPlusTwelve())
    c.write(join(test_dir, 'EightPlusTwelve.mol'))


def test_Dodecahedron():
    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = Cage([bb1, bb2], Dodecahedron())
    c.write(join(test_dir, 'Dodecahedron.mol'))


def test_multiconformer():
    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb1.update_from_mol(join(data_dir, 'amine2_conf2.mol'), 1)
    bb2 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    bb2.update_from_mol(join(data_dir, 'aldehyde3_conf2.mol'), 1)

    c = Cage([bb1, bb2],
             FourPlusSix(),
             bb_conformers=[0, 0])
    c.add_conformer([1, 0])
    c.add_conformer([0, 1])
    c.add_conformer([1, 1])

    c.write(join(test_dir, 'FourPlusSix_conf1.mol'), 0)
    c.write(join(test_dir, 'FourPlusSix_conf2.mol'), 1)
    c.write(join(test_dir, 'FourPlusSix_conf3.mol'), 2)
    c.write(join(test_dir, 'FourPlusSix_conf4.mol'), 3)
