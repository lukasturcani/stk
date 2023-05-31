import itertools as it

import pytest
import stk
from pytest_lazyfixture import lazy_fixture

# Fixtures need be visible for lazy_fixture() calls.
from .fixtures import *  # noqa


@pytest.fixture(
    params=(
        lazy_fixture("default_init"),
        lazy_fixture("init_from_file"),
        lazy_fixture("init"),
        lazy_fixture("init_from_molecule"),
        lazy_fixture("init_from_rdkit_mol"),
    ),
)
def case_data(request):
    """
    A :class:`.CaseData` instance.

    """

    return request.param


@pytest.fixture
def building_block(case_data):
    """
    A :class:`.BuildingBlock` instance.

    """

    return case_data.building_block


@pytest.fixture(
    params=(
        lambda molecule: stk.BromoFactory().get_functional_groups(molecule),
        lambda molecule: stk.PrimaryAminoFactory().get_functional_groups(
            molecule
        ),
        lambda molecule: it.chain(
            stk.PrimaryAminoFactory().get_functional_groups(molecule),
            stk.BromoFactory().get_functional_groups(molecule),
        ),
    )
)
def get_functional_groups(request):
    """
    Yield the functional groups of a `molecule`.

    Parameters
    ----------
    molecule : :class:`.Molecule`
        The molecule whose functional groups should be gotten.

    Yields
    ------
    :class:`.FunctionalGroup`
        A functional group of `molecule`.

    """

    return request.param
