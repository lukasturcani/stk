import itertools as it

import numpy as np
import pytest
import stk

from ..case_data import CaseData


@pytest.fixture
def functional_group2_2(functional_group2):
    """
    A :class:`.GenericFunctionalGroup` instance with 2 bonder atoms.

    """

    return functional_group2


@pytest.fixture
def two_two_reaction(
    functional_group2,
    functional_group2_2,
    bond_order,
    periodicity,
):
    """
    A :class:`.CaseData` instance.

    """

    return CaseData(
        reaction=stk.TwoTwoReaction(
            construction_state=get_construction_state(
                functional_group1=functional_group2_2,
                functional_group2=functional_group2,
            ),
            functional_group1=functional_group2_2,
            functional_group2=functional_group2,
            bond_order=bond_order,
            periodicity=periodicity,
        ),
        new_atoms=(),
        new_bonds=tuple(
            get_bonds(
                functional_group1=functional_group2_2,
                functional_group2=functional_group2,
                bond_order=bond_order,
                periodicity=periodicity,
            ),
        ),
        deleted_atoms=tuple(
            it.chain(
                functional_group2_2.get_deleters(),
                functional_group2.get_deleters(),
            )
        ),
        deleted_bonds=(),
    )


def get_bonds(
    functional_group1,
    functional_group2,
    bond_order,
    periodicity,
):
    b1, b2 = functional_group1.get_bonders()
    b3, b4 = functional_group2.get_bonders()
    yield stk.Bond(
        atom1=b1,
        atom2=b3,
        order=bond_order,
        periodicity=periodicity,
    )

    yield stk.Bond(
        atom1=b2,
        atom2=b4,
        order=bond_order,
        periodicity=periodicity,
    )


class MockConstructionState(stk.ConstructionState):
    def __init__(self, position_matrix):
        self.position_matrix = position_matrix

    def get_position_matrix(self):
        return self.position_matrix


def get_construction_state(functional_group1, functional_group2):
    size = max(
        it.chain(
            functional_group1.get_bonder_ids(),
            functional_group2.get_bonder_ids(),
        )
    )
    position_matrix = np.zeros((size + 1, 3))
    b1, b2 = functional_group1.get_bonder_ids()
    b3, b4 = functional_group2.get_bonder_ids()
    position_matrix[b1] = [100, 100, 100]
    position_matrix[b3] = [99, 99, 99]
    position_matrix[b2] = [10, 10, 10]
    position_matrix[b4] = [8, 8, 8]
    return MockConstructionState(position_matrix)
