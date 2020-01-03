import itertools as it
import pytest
import stk

from ._test_case import _TestCase


@pytest.fixture
def functional_group2(functional_group1):
    return functional_group1


@pytest.fixture
def one_one_reaction(
    functional_group1,
    functional_group2,
    bond_order,
    periodicity,
):
    return _TestCase(
        reaction=stk.OneOneReaction(
            functional_group1=functional_group1,
            functional_group2=functional_group2,
            bond_order=bond_order,
            periodicity=periodicity,
        ),
        new_atoms=(),
        new_bonds=(
            get_bond(
                functional_group1=functional_group1,
                functional_group2=functional_group2,
                bond_order=bond_order,
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
