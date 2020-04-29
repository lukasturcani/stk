import pytest
import itertools as it
import stk
from stk.molecular.reactions.reactions.reaction import ReactionResult

from ..case_data import CaseData
from .utilities import MockConstructionState, MockEdge


@pytest.fixture
def one_two_reaction(
    periodicity,
    functional_group1,
    functional_group2,
    bond_order,
    is_dative,
):
    bond_order_key = frozenset({
        type(functional_group1),
        type(functional_group2),
    })
    is_dative_key = frozenset({
        type(functional_group1),
        type(functional_group2),
    })
    edge = MockEdge(0, periodicity)
    return CaseData(
        factory=stk.GenericReactionFactory(
            bond_orders={
                bond_order_key: bond_order,
            },
            is_datives={
                is_dative_key: is_dative,
            },
        ),
        construction_state=MockConstructionState(
            edges=(edge, ),
            edge_functional_groups={
                0: (
                    functional_group1,
                    functional_group2,
                ),
            },
        ),
        edge_group=stk.EdgeGroup(
            edges=(edge, ),
        ),
        reaction_result=ReactionResult(
            new_atoms=(),
            new_bonds=get_new_bonds(
                functional_group1=functional_group1,
                functional_group2=functional_group2,
                order=bond_order,
                periodicity=periodicity,
                is_dative=is_dative,
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
    is_dative,
):
    fg1_bonder = next(functional_group1.get_bonders())
    fg2_bonders = functional_group2.get_bonders()
    yield stk.Bond(
        atom1=fg1_bonder,
        atom2=next(fg2_bonders),
        order=order,
        periodicity=periodicity,
        is_dative=is_dative,
    )
    yield stk.Bond(
        atom1=fg1_bonder,
        atom2=next(fg2_bonders),
        order=order,
        periodicity=periodicity,
        is_dative=is_dative,
    )
