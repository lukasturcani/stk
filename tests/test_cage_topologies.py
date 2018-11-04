import stk
import os
from os.path import join


test_dir = 'cage_topology_tests'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)
data_dir = os.path.join(os.getcwd(), 'data', 'cage_topologies')


# 3 + 4 topology tests.
def test_SixPlusEight():
    bb1 = stk.StructUnit3(join(data_dir, 'aldehyde3.mol'))
    bb2 = stk.StructUnit3(join(data_dir, 'amine4.mol'))
    c = stk.Cage([bb1, bb2], stk.SixPlusEight())
    c.write(join(test_dir, 'SixPlusEight.mol'))


# 2 + 4 topolgy tests.
def test_TwoPlusFour():
    bb1 = stk.StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = stk.StructUnit3(join(data_dir, 'amine4.mol'))
    c = stk.Cage([bb1, bb2], stk.TwoPlusFour())
    c.write(join(test_dir, 'TwoPlusFour.mol'))


def test_ThreePlusSix():
    bb1 = stk.StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = stk.StructUnit3(join(data_dir, 'amine4.mol'))
    c = stk.Cage([bb1, bb2], stk.ThreePlusSix())
    c.write(join(test_dir, 'ThreePlusSix.mol'))


def test_FourPlusEight():
    bb1 = stk.StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = stk.StructUnit3(join(data_dir, 'amine4.mol'))
    c = stk.Cage([bb1, bb2], stk.FourPlusEight())
    c.write(join(test_dir, 'FourPlusEight.mol'))


def test_FivePlusTen():
    bb1 = stk.StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = stk.StructUnit3(join(data_dir, 'amine4.mol'))
    c = stk.Cage([bb1, bb2], stk.FivePlusTen())
    c.write(join(test_dir, 'FivePlusTen.mol'))


def test_SixPlusTwelve():
    bb1 = stk.StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = stk.StructUnit3(join(data_dir, 'amine4.mol'))
    c = stk.Cage([bb1, bb2], stk.SixPlusTwelve())
    c.write(join(test_dir, 'SixPlusTwelve.mol'))


def test_EightPlusSixteen():
    bb1 = stk.StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = stk.StructUnit3(join(data_dir, 'amine4.mol'))
    c = stk.Cage([bb1, bb2], stk.EightPlusSixteen())
    c.write(join(test_dir, 'EightPlusSixteen.mol'))


def test_TenPlusTwenty():
    bb1 = stk.StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = stk.StructUnit3(join(data_dir, 'amine4.mol'))
    c = stk.Cage([bb1, bb2], stk.TenPlusTwenty())
    c.write(join(test_dir, 'TenPlusTwenty.mol'))


# 3 + 3 cage topologies.
def test_OnePlusOne():
    bb1 = stk.StructUnit3(join(data_dir, 'aldehyde3.mol'))
    bb2 = stk.StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = stk.Cage([bb1, bb2], stk.OnePlusOne())
    c.write(join(test_dir, 'OnePlusOne.mol'))


def test_TwoPlusTwo():
    bb1 = stk.StructUnit3(join(data_dir, 'aldehyde3.mol'))
    bb2 = stk.StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = stk.Cage([bb1, bb2], stk.TwoPlusTwo())
    c.write(join(test_dir, 'TwoPlusTwo.mol'))


def test_FourPlusFour():
    bb1 = stk.StructUnit3(join(data_dir, 'aldehyde3.mol'))
    bb2 = stk.StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = stk.Cage([bb1, bb2], stk.FourPlusFour())
    c.write(join(test_dir, 'FourPlusFour.mol'))


# 2 + 3 cage topologies.
def test_TwoPlusThree():
    bb1 = stk.StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = stk.StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = stk.Cage([bb1, bb2], stk.TwoPlusThree())
    c.write(join(test_dir, 'TwoPlusThree.mol'))


def test_FourPlusSix():
    bb1 = stk.StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = stk.StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = stk.Cage([bb1, bb2], stk.FourPlusSix())
    c.write(join(test_dir, 'FourPlusSix.mol'))


def test_multiFourPlusSix():
    bb1 = stk.StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = stk.StructUnit2(join(data_dir, 'amine2_1.mol'))
    bb3 = stk.StructUnit2(join(data_dir, 'amine2_2.mol'))

    bb4 = stk.StructUnit3(join(data_dir, 'aldehyde3.mol'))
    bb5 = stk.StructUnit3(join(data_dir, 'aldehyde3_1.mol'))
    bb6 = stk.StructUnit3(join(data_dir, 'aldehyde3_2.mol'))

    c1 = stk.Cage([bb1, bb2, bb3, bb4, bb5, bb6], stk.FourPlusSix())

    c2 = stk.Cage(
              [bb1, bb2, bb3, bb4, bb5, bb6],
              stk.FourPlusSix(
                  bb_assignments={
                      0: [0, 1],
                      1: [2, 3, 4],
                      2: [5],
                      3: [0],
                      4: [1, 2],
                      5: [3]
                   }
              )
    )

    c1.write(join(test_dir, 'multi_FourPlusSix_1.mol'))
    c2.write(join(test_dir, 'multi_FourPlusSix_2.mol'))


def test_FourPlusSix2():
    bb1 = stk.StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = stk.StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = stk.Cage([bb1, bb2], stk.FourPlusSix2())
    c.write(join(test_dir, 'FourPlusSix2.mol'))


def test_SixPlusNine():
    bb1 = stk.StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = stk.StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = stk.Cage([bb1, bb2], stk.SixPlusNine())
    c.write(join(test_dir, 'SixPlusNine.mol'))


def test_EightPlusTwelve():
    bb1 = stk.StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = stk.StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = stk.Cage([bb1, bb2], stk.EightPlusTwelve())
    c.write(join(test_dir, 'EightPlusTwelve.mol'))


def test_Dodecahedron():
    bb1 = stk.StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = stk.StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = stk.Cage([bb1, bb2], stk.Dodecahedron())
    c.write(join(test_dir, 'Dodecahedron.mol'))


def test_multiconformer():
    bb1 = stk.StructUnit2(join(data_dir, 'amine2.mol'))
    bb1.update_from_mol(join(data_dir, 'amine2_conf2.mol'), 1)
    bb2 = stk.StructUnit3(join(data_dir, 'aldehyde3.mol'))
    bb2.update_from_mol(join(data_dir, 'aldehyde3_conf2.mol'), 1)

    c = stk.Cage([bb1, bb2],
                 stk.FourPlusSix(),
                 bb_conformers=[0, 0])
    c.add_conformer([1, 0])
    c.add_conformer([0, 1])
    c.add_conformer([1, 1])

    c.write(join(test_dir, 'FourPlusSix_conf1.mol'), 0)
    c.write(join(test_dir, 'FourPlusSix_conf2.mol'), 1)
    c.write(join(test_dir, 'FourPlusSix_conf3.mol'), 2)
    c.write(join(test_dir, 'FourPlusSix_conf4.mol'), 3)
