import itertools as it

import pytest
import stk

from ..case_data import CaseData
from .utilities import MockConstructionState, MockEdge


@pytest.fixture
def one_one_reaction(
    periodicity,
    functional_group1,
    functional_group1_2,
    bond_order,
):
    bond_order_key = frozenset(
        {
            type(functional_group1),
            type(functional_group1_2),
        }
    )
    edge = MockEdge(0, periodicity)
    return CaseData(
        factory=stk.GenericReactionFactory(
            bond_orders={
                bond_order_key: bond_order,
            },
        ),
        construction_state=MockConstructionState(
            edges=(edge,),
            edge_functional_groups={
                0: (
                    functional_group1,
                    functional_group1_2,
                )
            },
        ),
        edge_group=stk.EdgeGroup(edges=(edge,)),
        reaction_result=stk.ReactionResult(
            new_atoms=(),
            new_bonds=get_new_bonds(
                functional_group1=functional_group1,
                functional_group2=functional_group1_2,
                order=bond_order,
                periodicity=periodicity,
            ),
            deleted_atoms=it.chain(
                functional_group1.get_deleters(),
                functional_group1_2.get_deleters(),
            ),
            deleted_bonds=(),
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
