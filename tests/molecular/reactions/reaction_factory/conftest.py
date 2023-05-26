import pytest
import stk
from pytest_lazyfixture import lazy_fixture

# Fixtures must be visible for lazy_fixture() calls.
from .fixtures import *  # noqa


@pytest.fixture(
    scope="session",
    params=(
        lambda: stk.Alcohol(
            oxygen=stk.O(0),
            hydrogen=stk.H(1),
            atom=stk.C(2),
            bonders=(stk.C(2),),
            deleters=(stk.O(0), stk.H(1)),
        ),
        lambda: stk.Aldehyde(
            carbon=stk.C(0),
            oxygen=stk.O(1),
            hydrogen=stk.H(2),
            atom=stk.C(3),
            bonders=(stk.C(0),),
            deleters=(stk.O(1),),
        ),
        lambda: stk.Alkene(
            carbon1=stk.C(0),
            atom1=stk.H(1),
            atom2=stk.H(2),
            carbon2=stk.C(3),
            atom3=stk.C(4),
            atom4=stk.H(5),
            bonders=(stk.C(4),),
            deleters=(stk.C(0), stk.H(1), stk.H(2)),
        ),
        lambda: stk.Alkyne(
            carbon1=stk.C(0),
            atom1=stk.H(1),
            carbon2=stk.C(2),
            atom2=stk.C(3),
            bonders=(stk.C(2),),
            deleters=(stk.C(0), stk.H(1)),
        ),
        lambda: stk.Amide(
            carbon=stk.C(0),
            oxygen=stk.O(1),
            nitrogen=stk.N(2),
            hydrogen1=stk.H(3),
            hydrogen2=stk.H(4),
            atom=stk.C(5),
            bonders=(stk.C(0),),
            deleters=(stk.N(2), stk.H(3), stk.H(4)),
        ),
        lambda: stk.BoronicAcid(
            boron=stk.B(0),
            oxygen1=stk.O(1),
            hydrogen1=stk.H(2),
            oxygen2=stk.O(3),
            hydrogen2=stk.H(4),
            atom=stk.C(5),
            bonders=(stk.B(0),),
            deleters=(stk.O(1), stk.H(2)),
        ),
        lambda: stk.Bromo(
            bromine=stk.Br(0),
            atom=stk.C(1),
            bonders=(stk.C(1),),
            deleters=(stk.Br(0),),
        ),
        lambda: stk.CarboxylicAcid(
            carbon=stk.C(0),
            oxygen1=stk.O(1),
            oxygen2=stk.O(2),
            hydrogen=stk.H(3),
            atom=stk.C(4),
            bonders=(stk.C(0),),
            deleters=(stk.O(2), stk.H(3)),
        ),
        lambda: stk.Dibromo(
            bromine1=stk.Br(0),
            atom1=stk.C(1),
            bromine2=stk.Br(2),
            atom2=stk.C(3),
            bonders=(stk.C(1),),
            deleters=(stk.Br(0),),
        ),
        lambda: stk.Difluoro(
            fluorine1=stk.F(0),
            atom1=stk.C(1),
            fluorine2=stk.F(2),
            atom2=stk.C(3),
            bonders=(stk.C(1),),
            deleters=(stk.F(0),),
        ),
        lambda: stk.Diol(
            atom1=stk.C(0),
            oxygen1=stk.O(1),
            hydrogen1=stk.H(2),
            atom2=stk.C(3),
            oxygen2=stk.O(4),
            hydrogen2=stk.H(5),
            bonders=(stk.C(0),),
            deleters=(stk.O(1), stk.H(2)),
        ),
        lambda: stk.Fluoro(
            fluorine=stk.F(1),
            atom=stk.C(0),
            bonders=(stk.C(0),),
            deleters=(stk.F(1),),
        ),
        lambda: stk.Iodo(
            iodine=stk.I(0),
            atom=stk.C(1),
            bonders=(stk.C(1),),
            deleters=(stk.I(0),),
        ),
        lambda: stk.PrimaryAmino(
            nitrogen=stk.N(0),
            hydrogen1=stk.H(1),
            hydrogen2=stk.H(2),
            atom=stk.C(3),
            bonders=(stk.C(3),),
            deleters=(stk.N(0), stk.H(1), stk.H(2)),
        ),
        lambda: stk.SecondaryAmino(
            nitrogen=stk.N(0),
            hydrogen=stk.H(1),
            atom1=stk.C(2),
            atom2=stk.C(3),
            bonders=(stk.N(0),),
            deleters=(stk.H(1),),
        ),
        lambda: stk.Thioacid(
            carbon=stk.C(0),
            oxygen=stk.O(1),
            sulfur=stk.S(2),
            hydrogen=stk.H(3),
            atom=stk.C(4),
            bonders=(stk.C(0),),
            deleters=(stk.S(2), stk.H(3)),
        ),
        lambda: stk.Thiol(
            sulfur=stk.S(0),
            hydrogen=stk.H(1),
            atom=stk.C(2),
            bonders=(stk.C(2),),
            deleters=(stk.S(0), stk.H(1)),
        ),
        lambda: stk.SingleAtom(atom=stk.Fe(0)),
    ),
)
def functional_group1(request) -> stk.GenericFunctionalGroup:
    """
    A :class:`.GenericFunctionalGroup` with 1 bonder atom.

    """

    return request.param()


@pytest.fixture
def functional_group1_2(functional_group1):
    """
    A :class:`.GenericFunctionalGroup` with 1 bonder atom.

    """

    return functional_group1


@pytest.fixture(
    params=(
        lambda: stk.Dibromo(
            bromine1=stk.Br(0),
            atom1=stk.C(1),
            bromine2=stk.Br(2),
            atom2=stk.C(3),
            bonders=(stk.C(1), stk.C(2)),
            deleters=(stk.Br(0), stk.Br(2)),
        ),
        lambda: stk.Difluoro(
            fluorine1=stk.F(0),
            atom1=stk.C(1),
            fluorine2=stk.F(2),
            atom2=stk.C(3),
            bonders=(stk.C(1), stk.C(3)),
            deleters=(stk.F(0), stk.F(2)),
        ),
        lambda: stk.Diol(
            atom1=stk.C(0),
            oxygen1=stk.O(1),
            hydrogen1=stk.H(2),
            atom2=stk.C(3),
            oxygen2=stk.O(4),
            hydrogen2=stk.H(5),
            bonders=(stk.C(0), stk.C(3)),
            deleters=(stk.O(1), stk.H(2), stk.O(4), stk.H(5)),
        ),
    ),
)
def functional_group2(request):
    """
    A :class:`.GenericFunctionalGroup` with 2 bonder atoms.

    """

    return request.param()


@pytest.fixture(
    params=(
        (0, 0, 0),
        (-1, 0, 1),
        (-1, -1, -2),
        (1, 2, 3),
    ),
)
def periodicity(request):
    """
    The periodicity of a bond created by a :class:`.Reaction`.

    """

    return request.param


@pytest.fixture(
    params=(1, 2, 9),
)
def bond_order(request):
    """
    The bond order of a bond created by a :class:`.Reaction`.

    """

    return request.param


@pytest.fixture(
    params=(
        lazy_fixture("one_one_reaction"),
        lazy_fixture("one_two_reaction"),
        lazy_fixture("two_two_reaction"),
        lazy_fixture("dative_reaction"),
    ),
)
def case_data(request):
    """
    A :class:`.CaseData` instance.

    """

    return request.param
