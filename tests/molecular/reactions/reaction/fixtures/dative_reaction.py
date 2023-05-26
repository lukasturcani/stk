import itertools as it

import pytest
import stk

from ..case_data import CaseData


@pytest.fixture
def dative_reaction(
    functional_group1,
    functional_group1_2,
    bond_order,
    periodicity,
):
    """
    A :class:`.CaseData` instance.

    """

    return CaseData(
        reaction=stk.DativeReaction(
            reaction=stk.OneOneReaction(
                functional_group1=functional_group1,
                functional_group2=functional_group1_2,
                bond_order=bond_order,
                periodicity=periodicity,
            ),
        ),
        new_atoms=(),
        new_bonds=(
            get_bond(
                functional_group1=functional_group1,
                functional_group2=functional_group1_2,
                bond_order=bond_order,
                periodicity=periodicity,
            ),
        ),
        deleted_atoms=tuple(
            it.chain(
                functional_group1.get_deleters(),
                functional_group1_2.get_deleters(),
            )
        ),
        deleted_bonds=(),
    )


def get_bond(
    functional_group1,
    functional_group2,
    bond_order,
    periodicity,
):
    (bonder1,) = functional_group1.get_bonders()
    (bonder2,) = functional_group2.get_bonders()

    if is_metal(bonder1):
        return stk.Bond(
            atom1=bonder2,
            atom2=bonder1,
            order=bond_order,
            periodicity=periodicity,
        )
    else:
        return stk.Bond(
            atom1=bonder1,
            atom2=bonder2,
            order=bond_order,
            periodicity=periodicity,
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
