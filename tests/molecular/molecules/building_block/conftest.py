import itertools as it

import pytest
import stk
from pytest_lazyfixture import lazy_fixture
from collections import abc

# Fixtures need be visible for lazy_fixture() calls.
from .fixtures import *  # noqa

from .case_data import CaseData


@pytest.fixture(
    params=(
        lazy_fixture("default_init"),
        lazy_fixture("init_from_file"),
        lazy_fixture("init"),
        lazy_fixture("init_from_molecule"),
        lazy_fixture("init_from_rdkit_mol"),
    ),
)
def case_data(request) -> CaseData:
    """
    A :class:`.CaseData` instance.

    """

    return request.param


@pytest.fixture
def building_block(case_data: CaseData) -> stk.BuildingBlock:
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
def get_functional_groups(
    request: pytest.FixtureRequest,
) -> abc.Iterable[stk.FunctionalGroup]:
    """
    Get the functional groups of a `molecule`.

    Parameters:
        molecule:
            The molecule whose functional groups should be gotten.

    Returns:
        A functional group of `molecule`.

    """

    return request.param
