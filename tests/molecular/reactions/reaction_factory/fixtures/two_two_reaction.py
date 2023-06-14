import itertools as it

import numpy as np
import pytest
import stk

from ..case_data import CaseData
from .utilities import MockConstructionState, MockEdge


@pytest.fixture
def functional_group2_2(functional_group2):
    """
    A :class:`.GenericFunctionalGroup` with 2 bonder atoms.

    """

    return functional_group2.clone()


@pytest.fixture
def two_two_reaction(
    periodicity,
    functional_group2,
    functional_group2_2,
    bond_order,
):
    bond_order_key = frozenset(
        {
            type(functional_group2),
            type(functional_group2_2),
        }
    )
    position_matrix = get_position_matrix(
        functional_group1=functional_group2,
        functional_group2=functional_group2_2,
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
                    functional_group2,
                    functional_group2_2,
                ),
            },
            position_matrix=position_matrix,
        ),
        edge_group=stk.EdgeGroup(
            edges=(edge,),
        ),
        reaction_result=stk.ReactionResult(
            new_atoms=(),
            new_bonds=get_new_bonds(
                functional_group1=functional_group2,
                functional_group2=functional_group2_2,
                order=bond_order,
                periodicity=periodicity,
            ),
            deleted_atoms=it.chain(
                functional_group2.get_deleters(),
                functional_group2_2.get_deleters(),
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
    bonder1, bonder2 = functional_group1.get_bonders()
    bonder3, bonder4 = functional_group2.get_bonders()
    yield stk.Bond(
        atom1=bonder1,
        atom2=bonder3,
        order=order,
        periodicity=periodicity,
    )
    yield stk.Bond(
        atom1=bonder2,
        atom2=bonder4,
        order=order,
        periodicity=periodicity,
    )


def get_position_matrix(functional_group1, functional_group2):
    size = (
        max(
            it.chain(
                functional_group1.get_bonder_ids(),
                functional_group2.get_bonder_ids(),
            )
        )
        + 1
    )
    position_matrix = np.zeros((size, 3))
    bonder1, bonder2 = functional_group1.get_bonder_ids()
    bonder3, bonder4 = functional_group2.get_bonder_ids()

    # Place the first bonder of each functional group close.
    position_matrix[bonder1] = [100, 100, 100]
    position_matrix[bonder3] = [99, 100, 100]
    return position_matrix
