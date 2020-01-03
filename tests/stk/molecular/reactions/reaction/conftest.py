import pytest
import stk

from .fixtures import *


@pytest.fixture(
    params=(
        stk.Alcohol(
            oxygen=stk.O(1),
            hydrogen=stk.H(2),
            atom=stk.C(2),
            bonders=(stk.C(2), ),
            deleters=(stk.O(1), stk.H(2)),
        ),

        stk.Aldehyde(
            carbon=stk.C(0),
            oxygen=stk.O(1),
            hydrogen=stk.H(2),
            atom=stk.C(3),
            bonders=(stk.C(0), ),
            deleters=(stk.O(1), ),
        ),

        stk.Alkene(
            carbon1=stk.C(1),
            atom1=stk.C(2),
            atom2=stk.C(3),
            carbon2=stk.C(4),
            atom3=stk.C(5),
            atom4=stk.C(6),
            bonders=(stk.C(1), ),
            deleters=(stk.C(5), stk.C(6)),
        ),

        stk.Alkyne(
            carbon1=stk.C(0),
            atom1=stk.C(1),
            carbon2=stk.C(2),
            atom2=stk.H(3),
            bonders=(stk.C(0), ),
            deleters=(stk.C(2), stk.H(3)),
        ),

        stk.Amide(
            carbon=stk.C(0),
            oxygen=stk.O(1),
            nitrogen=stk.N(2),
            hydrogen1=stk.H(3),
            hydrogen2=stk.H(4),
            atom=stk.C(5),
            bonders=(stk.C(0), ),
            deleters=(stk.N(2), stk.H(3), stk.H(4)),
        ),

        stk.BoronicAcid(
            boron=stk.B(0),
            oxygen1=stk.O(1),
            hydrogen1=stk.H(2),
            oxygen2=stk.O(3),
            hydrogen2=stk.H(4),
            atom=stk.C(5),
            bonders=(stk.B(0), ),
            deleters=(stk.O(1), stk.H(2), stk.O(3), stk.H(4)),
        ),

        stk.Bromo(
            bromine=stk.Br(0),
            atom=stk.C(1),
            bonders=(stk.C(1), ),
            deleters=(stk.Br(0), ),
        ),

        stk.CarboxylicAcid(
            carbon=stk.C(0),
            oxygen1=stk.O(1),
            oxygen2=stk.O(5),
            hydrogen=stk.H(6),
            atom=stk.C(8),
            bonders=(stk.C(0), ),
            deleters=(stk.O(5), stk.H(6)),
        ),

        stk.Fluoro(
            fluorine=stk.F(12),
            atom=stk.C(2),
            bonders=(stk.C(2), ),
            deleters=(stk.F(12), ),
        ),

        stk.Iodo(
            iodine=stk.I(12),
            atom=stk.C(0),
            bonders=(stk.C(0), ),
            deleters=(stk.I(12), ),
        ),

        stk.PrimaryAmino(
            nitrogen=stk.N(0),
            hydrogen1=stk.H(2),
            hydrogen2=stk.H(3),
            atom=stk.C(4),
            bonders=(stk.C(4), ),
            deleters=(stk.N(0), stk.H(2), stk.H(3)),
        ),

        stk.SecondaryAmino(
            nitrogen=stk.N(12),
            hydrogen=stk.H(21),
            atom1=stk.C(0),
            atom2=stk.C(12),
            bonders=(stk.N(12), ),
            deleters=(stk.H(21), ),
        ),

        stk.Thioacid(
            carbon=stk.C(0),
            oxygen=stk.O(12),
            sulfur=stk.S(21),
            hydrogen=stk.H(2),
            atom=stk.C(1),
        ),

        stk.Thiol(
            sulfur=stk.S(12),
            hydrogen=stk.H(21),
            atom=stk.C(32),
            bonders=(stk.C(32), ),
            deleters=(stk.S(12), stk.H(21)),
        ),

    ),
)
def functional_group1(request):
    return request.param


@pytest.fixture(
    params=(

        stk.BoronicAcid(
            boron=stk.B(0),
            oxygen1=stk.O(1),
            hydrogen1=stk.H(2),
            oxygen2=stk.O(3),
            hydrogen2=stk.H(4),
            atom=stk.C(5),
            bonders=(stk.O(1), stk.O(3)),
            deleters=(stk.H(2), stk.H(4)),
        ),

        stk.Dibromo(
            bromine1=stk.Br(0),
            atom1=stk.C(1),
            bromine2=stk.Br(32),
            atom2=stk.C(31),
            bonders=(stk.C(1), stk.C(31)),
            deleters=(stk.Br(0), stk.Br(32)),
        ),

        stk.Difluoro(
            fluorine1=stk.F(0),
            atom1=stk.C(1),
            fluorine2=stk.F(3),
            atom2=stk.C(4),
            bonders=(stk.C(1), stk.C(4)),
            deleters=(stk.F(0), stk.F(3)),
        ),

        stk.Diol(
            atom1=stk.C(1),
            oxygen1=stk.O(2),
            hydrogen1=stk.H(0),
            atom2=stk.C(32),
            oxygen2=stk.O(33),
            hydrogen2=stk.H(43),
            bonders=(stk.C(1), stk.C(32)),
            deleters=(stk.O(2), stk.H(0), stk.O(33), stk.H(43)),
        ),

    ),
)
def functional_group2(request):
    return request.param


@pytest.fixture(
    params=(
        pytest.lazy_fixture('one_one_reaction'),
        pytest.lazy_fixture('one_two_reaction'),
        pytest.lazy_fixture('two_two_reaction'),
        pytest.lazy_fixture('ring_amine_reaction'),
    ),
)
def test_case(request):
    return request.param


@pytest.fixture(
    params=(1, 2, 3),
)
def bond_order(request):
    return request.param


@pytest.fixture(
    params=(
        (0, 0, 0),
        (1, 0, 0),
        (0, -1, 0),
    ),
)
def periodicity(request):
    return request.param
