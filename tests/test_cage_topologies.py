from ..molecular.topologies.cage import *
from ..molecular import StructUnit3, StructUnit2, Cage
import os
from os.path import join

test_dir = 'cage_topology_tests'
os.mkdir(test_dir)
data_dir = os.path.join(os.path.getcwd(), 'data', 'cage_topologies')

# 3 + 4 topology tests.
def test_SixPlusEight():
    if test_dir not in os.getcwd():
        os.chdir(test_dir)

    bb1 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    bb2 = StructUnit3(join(data_dir, 'amine4.mol'))
    c = Cage([bb1, bb2], SixPlusEight())

# 2 + 4 topolgy tests.
def test_TwoPlusFour():
    if test_dir not in os.getcwd():
        os.chdir(test_dir)

    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = Cage([bb1, bb2], SixPlusEight())


def test_ThreePlusSix():
    if test_dir not in os.getcwd():
        os.chdir(test_dir)

    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = Cage([bb1, bb2], SixPlusEight())


def test_FourPlusEight():
    if test_dir not in os.getcwd():
        os.chdir(test_dir)

    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = Cage([bb1, bb2], SixPlusEight())


def test_SixPlusTwelve():
    if test_dir not in os.getcwd():
        os.chdir(test_dir)

    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = Cage([bb1, bb2], SixPlusEight())


def test_TenPlusTwenty():
    if test_dir not in os.getcwd():
        os.chdir(test_dir)

    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = Cage([bb1, bb2], SixPlusEight())

# 3 + 3 cage topologies.
def test_OnePlusOne():
    if test_dir not in os.getcwd():
        os.chdir(test_dir)

    bb1 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    bb2 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = Cage([bb1, bb2], SixPlusEight())


def test_TwoPlusTwo():
    if test_dir not in os.getcwd():
        os.chdir(test_dir)

    bb1 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    bb2 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = Cage([bb1, bb2], SixPlusEight())


def test_FourPlusFour():
    if test_dir not in os.getcwd():
        os.chdir(test_dir)

    bb1 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    bb2 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = Cage([bb1, bb2], SixPlusEight())


# 2 + 3 cage topologies.
def test_TwoPlusThree():
    if test_dir not in os.getcwd():
        os.chdir(test_dir)

    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = Cage([bb1, bb2], SixPlusEight())


def test_FourPlusSix():
    if test_dir not in os.getcwd():
        os.chdir(test_dir)

    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = Cage([bb1, bb2], SixPlusEight())


def test_FourPlusSix2():
    if test_dir not in os.getcwd():
        os.chdir(test_dir)

    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = Cage([bb1, bb2], SixPlusEight())


def test_SixPlusNine():
    if test_dir not in os.getcwd():
        os.chdir(test_dir)

    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = Cage([bb1, bb2], SixPlusEight())


def test_EightPlusTwelve():
    if test_dir not in os.getcwd():
        os.chdir(test_dir)

    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = Cage([bb1, bb2], SixPlusEight())


def test_Dodecahedron():
    if test_dir not in os.getcwd():
        os.chdir(test_dir)

    bb1 = StructUnit2(join(data_dir, 'amine2.mol'))
    bb2 = StructUnit3(join(data_dir, 'aldehyde3.mol'))
    c = Cage([bb1, bb2], SixPlusEight())
