import pytest
import itertools as it
import stk

from ._test_case import _TestCase
from .utilities import MockEdge


@pytest.fixture(
    params=(
        stk.Alcohol(),
        stk.Aldehyde(),
        stk.Alkene(),
        stk.Alkyne(),
        stk.Amide(),
        stk.BoronicAcid(),
        stk.Bromo(),
        stk.CarboxylicAcid(),
        stk.Dibromo(),
        stk.Difluoro(),
        stk.Diol(),
        stk.Fluoro(),
        stk.Iodo(),
        stk.PrimaryAmino(),
        stk.SecondaryAmino(),
        stk.Thioacid(),
        stk.Thiol(),
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
        reaction_result=stk.ReactionResult(
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
