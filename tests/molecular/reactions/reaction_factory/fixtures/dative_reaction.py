import itertools as it

import pytest
import stk

from ..case_data import CaseData
from .utilities import MockConstructionState, MockEdge


@pytest.fixture
def dative_reaction(
    functional_group1,
    functional_group1_2,
    periodicity,
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
        factory=stk.DativeReactionFactory(
            stk.GenericReactionFactory(
                bond_orders={
                    bond_order_key: bond_order,
                },
            )
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


def is_metal(atom):
    metal_atomic_numbers = set(
        it.chain(
            range(21, 31),
            range(39, 49),
            range(72, 81),
        )
    )

    return atom.get_atomic_number() in metal_atomic_numbers


def get_new_bonds(
    functional_group1,
    functional_group2,
    order,
    periodicity,
):
    (bonder1,) = functional_group1.get_bonders()
    (bonder2,) = functional_group2.get_bonders()

    if is_metal(bonder1):
        yield stk.Bond(
            atom1=bonder2,
            atom2=bonder1,
            order=order,
            periodicity=periodicity,
        )
    else:
        yield stk.Bond(
            atom1=bonder1,
            atom2=bonder2,
            order=order,
            periodicity=periodicity,
        )
