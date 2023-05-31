import pytest
import stk
from pytest_lazyfixture import lazy_fixture

# Fixtures need to visible for lazy_fixture() calls.
from .fixtures import *  # noqa


@pytest.fixture(
    scope="session",
    params=(
        lambda: stk.Alcohol(
            oxygen=stk.O(1),
            hydrogen=stk.H(2),
            atom=stk.C(2),
            bonders=(stk.C(2),),
            deleters=(stk.O(1), stk.H(2)),
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
            carbon1=stk.C(1),
            atom1=stk.C(2),
            atom2=stk.C(3),
            carbon2=stk.C(4),
            atom3=stk.C(5),
            atom4=stk.C(6),
            bonders=(stk.C(1),),
            deleters=(stk.C(5), stk.C(6)),
        ),
        lambda: stk.Alkyne(
            carbon1=stk.C(0),
            atom1=stk.C(1),
            carbon2=stk.C(2),
            atom2=stk.H(3),
            bonders=(stk.C(0),),
            deleters=(stk.C(2), stk.H(3)),
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
            deleters=(stk.O(1), stk.H(2), stk.O(3), stk.H(4)),
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
            oxygen2=stk.O(5),
            hydrogen=stk.H(6),
            atom=stk.C(8),
            bonders=(stk.C(0),),
            deleters=(stk.O(5), stk.H(6)),
        ),
        lambda: stk.Fluoro(
            fluorine=stk.F(12),
            atom=stk.C(2),
            bonders=(stk.C(2),),
            deleters=(stk.F(12),),
        ),
        lambda: stk.Iodo(
            iodine=stk.I(12),
            atom=stk.C(0),
            bonders=(stk.C(0),),
            deleters=(stk.I(12),),
        ),
        lambda: stk.PrimaryAmino(
            nitrogen=stk.N(0),
            hydrogen1=stk.H(2),
            hydrogen2=stk.H(3),
            atom=stk.C(4),
            bonders=(stk.C(4),),
            deleters=(stk.N(0), stk.H(2), stk.H(3)),
        ),
        lambda: stk.SecondaryAmino(
            nitrogen=stk.N(12),
            hydrogen=stk.H(21),
            atom1=stk.C(0),
            atom2=stk.C(12),
            bonders=(stk.N(12),),
            deleters=(stk.H(21),),
        ),
        lambda: stk.Thioacid(
            carbon=stk.C(0),
            oxygen=stk.O(12),
            sulfur=stk.S(21),
            hydrogen=stk.H(2),
            atom=stk.C(1),
            bonders=(stk.C(0),),
            deleters=(stk.S(21), stk.H(2)),
        ),
        lambda: stk.Thiol(
            sulfur=stk.S(12),
            hydrogen=stk.H(21),
            atom=stk.C(32),
            bonders=(stk.C(32),),
            deleters=(stk.S(12), stk.H(21)),
        ),
        lambda: stk.SingleAtom(atom=stk.Fe(0)),
    ),
)
def functional_group1(request) -> stk.GenericFunctionalGroup:
    """
    A :class:`.GenericFunctionalGroup` instance with 1 bonder atom.

    """

    return request.param()


@pytest.fixture
def functional_group1_2(functional_group1):
    """
    A :class:`.GenericFunctionalGroup` instance with 1 bonder atom.

    """

    return functional_group1


@pytest.fixture(
    params=(
        lambda: stk.BoronicAcid(
            boron=stk.B(0),
            oxygen1=stk.O(1),
            hydrogen1=stk.H(2),
            oxygen2=stk.O(3),
            hydrogen2=stk.H(4),
            atom=stk.C(5),
            bonders=(stk.O(1), stk.O(3)),
            deleters=(stk.H(2), stk.H(4)),
        ),
        lambda: stk.Dibromo(
            bromine1=stk.Br(0),
            atom1=stk.C(1),
            bromine2=stk.Br(32),
            atom2=stk.C(31),
            bonders=(stk.C(1), stk.C(31)),
            deleters=(stk.Br(0), stk.Br(32)),
        ),
        lambda: stk.Difluoro(
            fluorine1=stk.F(0),
            atom1=stk.C(1),
            fluorine2=stk.F(3),
            atom2=stk.C(4),
            bonders=(stk.C(1), stk.C(4)),
            deleters=(stk.F(0), stk.F(3)),
        ),
        lambda: stk.Diol(
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
def functional_group2(request) -> stk.GenericFunctionalGroup:
    """
    A :class:`.GenericFunctionalGroup` instance with 2 bonder atoms.

    """

    return request.param()


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


@pytest.fixture(
    params=(1, 3, 9),
)
def bond_order(request):
    """
    The bond order of a bond made by a reaction.

    """

    return request.param


@pytest.fixture(
    params=(
        (0, 0, 0),
        (0, -1, 0),
    ),
)
def periodicity(request):
    """
    The periodicity of a bond made by a reaction.

    """

    return request.param
