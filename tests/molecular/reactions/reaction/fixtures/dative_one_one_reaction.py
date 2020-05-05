import itertools as it
import pytest
import stk

from ..case_data import CaseData


@pytest.fixture
def dative_one_one_reaction(
    dative_functional_groups,
    dative_bond_order,
    periodicity,
):
    """
    A :class:`.CaseData` instance.

    """

    functional_group1, functional_group2 = dative_functional_groups

    return CaseData(
        reaction=stk.OneOneReaction(
            functional_group1=functional_group1,
            functional_group2=functional_group2,
            bond_order=dative_bond_order,
            periodicity=periodicity,
        ),
        new_atoms=(),
        new_bonds=(
            get_bond(
                functional_group1=functional_group1,
                functional_group2=functional_group2,
                bond_order=dative_bond_order,
                periodicity=periodicity,
            ),
        ),
        deleted_atoms=tuple(it.chain(
            functional_group1.get_deleters(),
            functional_group2.get_deleters(),
        )),
    )


def get_bond(
    functional_group1,
    functional_group2,
    bond_order,
    periodicity,
):
    return stk.Bond(
        atom1=next(functional_group1.get_bonders()),
        atom2=next(functional_group2.get_bonders()),
        order=bond_order,
        periodicity=periodicity,
    )


def is_metal_atom(atom):
    # Metal atomic numbers.
    metal_atomic_numbers = set(it.chain(
        list(range(21, 31)),
        list(range(39, 49)),
        list(range(72, 81))
    ))

    return atom.get_atomic_number() in metal_atomic_numbers
