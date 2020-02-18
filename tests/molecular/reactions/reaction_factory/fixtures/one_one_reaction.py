import pytest
import itertools as it
import stk
from stk.molecular.reactions.reactions.reaction import ReactionResult

from ._test_case import _TestCase
from .utilities import MockEdge


@pytest.fixture(
    params=(
        stk.Alcohol(
            oxygen=stk.O(0),
            hydrogen=stk.H(1),
            atom=stk.C(2),
            bonders=(stk.C(2), ),
            deleters=(stk.O(0), stk.H(1)),
        ),
        stk.Aldehyde(
            carbon=stk.C(0),
            oxygen=stk.O(1),
            hydrogen=stk.H(2),
            atom=stk.C(3),
            bonders=(stk.C(0),),
            deleters=(stk.O(1),),
        ),
        stk.Alkene(
            carbon1=stk.C(0),
            atom1=stk.H(1),
            atom2=stk.H(2),
            carbon2=stk.C(3),
            atom3=stk.C(4),
            atom4=stk.H(5),
            bonders=(stk.C(4),),
            deleters=(stk.C(0), stk.H(1), stk.H(2)),
        ),
        stk.Alkyne(
            carbon1=stk.C(0),
            atom1=stk.H(1),
            carbon2=stk.C(2),
            atom2=stk.C(3),
            bonders=(stk.C(2),),
            deleters=(stk.C(0), stk.H(1)),
        ),
        stk.Amide(
            carbon=stk.C(0),
            oxygen=stk.O(1),
            nitrogen=stk.N(2),
            hydrogen1=stk.H(3),
            hydrogen2=stk.H(4),
            atom=stk.C(5),
            bonders=(stk.C(0),),
            deleters=(stk.N(2), stk.H(3), stk.H(4)),
        ),
        stk.BoronicAcid(
            boron=stk.B(0),
            oxygen1=stk.O(1),
            hydrogen1=stk.H(2),
            oxygen2=stk.O(3),
            hydrogen2=stk.H(4),
            atom=stk.C(5),
            bonders=(stk.B(0),),
            deleters=(stk.O(1), stk.H(2)),
        ),
        stk.Bromo(
            bromine=stk.Br(0),
            atom=stk.C(1),
            bonders=(stk.C(1),),
            deleters=(stk.Br(0),),
        ),
        stk.CarboxylicAcid(
            carbon=stk.C(0),
            oxygen1=stk.O(1),
            oxygen2=stk.O(2),
            hydrogen=stk.H(3),
            atom=stk.C(4),
            bonders=(stk.C(0), ),
            deleters=(stk.O(2), stk.H(3)),
        ),
        stk.Dibromo(
            bromine1=stk.Br(0),
            atom1=stk.C(1),
            bromine2=stk.Br(2),
            atom2=stk.C(2),
            bonders=(stk.C(1), ),
            deleters=(stk.Br(0), ),
        ),
        stk.Difluoro(
            fluorine1=stk.F(0),
            atom1=stk.C(1),
            fluorine2=stk.F(2),
            atom2=stk.C(3),
            bonders=(stk.C(1), ),
            deleters=(stk.F(0), ),
        ),
        stk.Diol(
            atom1=stk.C(0),
            oxygen1=stk.O(1),
            hydrogen1=stk.H(2),
            atom2=stk.C(3),
            oxygen2=stk.O(4),
            hydrogen2=stk.H(5),
            bonders=(stk.C(0), ),
            deleters=(stk.O(1), stk.H(2)),
        ),
        stk.Fluoro(
            fluorine=stk.F(1),
            atom=stk.C(0),
            bonders=(stk.C(0), ),
            deleters=(stk.F(1), ),
        ),
        stk.Iodo(
            iodine=stk.I(0),
            atom=stk.C(1),
            bonders=(stk.C(1), ),
            deleters=(stk.I(0), ),
        ),
        stk.PrimaryAmino(
            nitrogen=stk.N(0),
            hydrogen1=stk.H(1),
            hydrogen2=stk.H(2),
            atom=stk.C(3),
            bonders=(stk.C(3), ),
            deleters=(stk.N(0), stk.H(1), stk.H(2)),
        ),
        stk.SecondaryAmino(
            nitrogen=stk.N(0),
            hydrogen=stk.H(1),
            atom1=stk.C(2),
            atom2=stk.C(3),
            bonders=(stk.N(0), ),
            deleters=(stk.H(1), ),
        ),
        stk.Thioacid(
            carbon=stk.C(0),
            oxygen=stk.O(1),
            sulfur=stk.S(2),
            hydrogen=stk.H(3),
            atom=stk.C(4),
            bonders=(stk.C(0), ),
            deleters=(stk.S(2), stk.H(3)),
        ),
        stk.Thiol(
            sulfur=stk.S(0),
            hydrogen=stk.H(1),
            atom=stk.C(2),
            bonders=(stk.C(2), ),
            deleters=(stk.S(0), stk.H(1)),
        ),
    ),
)
def functional_group(request):
    return request.param


@pytest.fixture
def functional_group1(functional_group):
    return functional_group.clone()


@pytest.fixture
def functional_group2(functional_group):
    return functional_group.clone()


@pytest.fixture(
    params=(
        (0, 0, 0),
        (-1, 0, 1),
        (-1, -1, -2),
        (1, 2, 3),
    ),
)
def periodicity(request):
    return request.param


@pytest.fixture(
    params=(1, 2),
)
def bond_order(request):
    return request.param


@pytest.fixture
def one_one_reaction(
    periodicity,
    functional_group1,
    functional_group2,
    bond_order,
):
    bond_order_key = frozenset({
        type(functional_group1),
        type(functional_group2),
    })
    return _TestCase(
        factory=stk.GenericReactionFactory(
            bond_orders={
                bond_order_key: bond_order,
            },
        ),
        construction_state=None,
        edge=MockEdge(periodicity),
        functional_groups=(
            functional_group1,
            functional_group2,
        ),
        reaction_result=ReactionResult(
            new_atoms=(),
            new_bonds=get_new_bonds(
                functional_group1=functional_group1,
                functional_group2=functional_group2,
                order=bond_order,
                periodicity=periodicity,
            ),
            deleted_atoms=it.chain(
                functional_group1.get_deleters(),
                functional_group2.get_deleters(),
            ),
        ),
    )


def get_new_bonds(
    functional_group1,
    functional_group2,
    order,
    periodicity,
):
    yield stk.Bond(
        atom1=next(functional_group1.get_bonders()),
        atom2=next(functional_group2.get_bonders()),
        order=order,
        periodicity=periodicity,
    )
