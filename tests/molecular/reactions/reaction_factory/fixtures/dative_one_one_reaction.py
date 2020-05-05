import pytest
import itertools as it
import stk
from stk.molecular.reactions.reactions.reaction import ReactionResult

from ..case_data import CaseData
from .utilities import MockConstructionState, MockEdge


@pytest.fixture
def dative_one_one_reaction(
    periodicity,
    dative_functional_groups,
    dative_bond_order,
):
    functional_group1, functional_group2 = dative_functional_groups
    bond_order_key = frozenset({
        type(functional_group1),
        type(functional_group2),
    })
    edge = MockEdge(0, periodicity)
    return CaseData(
        factory=stk.DativeReactionFactory(
            bond_orders={
                bond_order_key: dative_bond_order,
            },
        ),
        construction_state=MockConstructionState(
            edges=(edge, ),
            edge_functional_groups={
                0: (
                    functional_group1,
                    functional_group2,
                )
            }
        ),
        edge_group=stk.EdgeGroup(
            edges=(edge, )
        ),
        reaction_result=ReactionResult(
            new_atoms=(),
            new_bonds=get_new_bonds(
                functional_group1=functional_group1,
                functional_group2=functional_group2,
                order=dative_bond_order,
                periodicity=periodicity,
            ),
            deleted_atoms=it.chain(
                functional_group1.get_deleters(),
                functional_group2.get_deleters(),
            ),
        ),
    )


def is_metal_atom(atom):
    # Metal atomic numbers.
    metal_atomic_numbers = set(it.chain(
        list(range(21, 31)),
        list(range(39, 49)),
        list(range(72, 81))
    ))

    return atom.get_atomic_number() in metal_atomic_numbers


def get_new_bonds(
    functional_group1,
    functional_group2,
    order,
    periodicity,
):

    bondera = next(functional_group1.get_bonders())
    bonderb = next(functional_group2.get_bonders())

    if is_metal_atom(bondera):
        bonder2 = bondera
        bonder1 = bonderb
    elif is_metal_atom(bonderb):
        bonder2 = bonderb
        bonder1 = bondera

    yield stk.Bond(
        atom1=bonder1,
        atom2=bonder2,
        order=order,
        periodicity=periodicity,
    )
