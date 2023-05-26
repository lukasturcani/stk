import itertools as it

import pytest
import stk

from ..case_data import CaseData


@pytest.fixture
def one_two_reaction(
    functional_group1,
    functional_group2,
    bond_order,
    periodicity,
):
    """
    A :class:`.CaseData` instance.

    """

    return CaseData(
        reaction=stk.OneTwoReaction(
            functional_group1=functional_group1,
            functional_group2=functional_group2,
            bond_order=bond_order,
            periodicity=periodicity,
        ),
        new_atoms=(),
        new_bonds=tuple(
            get_bonds(
                functional_group1=functional_group1,
                functional_group2=functional_group2,
                bond_order=bond_order,
                periodicity=periodicity,
            ),
        ),
        deleted_atoms=tuple(
            it.chain(
                functional_group1.get_deleters(),
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
    bonders1 = functional_group1.get_bonders()
    bonders2 = functional_group2.get_bonders()
    for bonder1, bonder2 in it.product(bonders1, bonders2):
        yield stk.Bond(
            atom1=bonder1,
            atom2=bonder2,
            order=bond_order,
            periodicity=periodicity,
        )
