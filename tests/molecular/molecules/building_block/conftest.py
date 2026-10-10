import itertools as it
from collections import abc

import pytest
from pytest_lazy_fixtures import lf

import stk

from .case_data import CaseData

# Fixtures need be visible for lazy_fixture() calls.
from .fixtures import *


@pytest.fixture(
    params=(
        lf("default_init"),
        lf("init_from_file"),
        lf("init"),
        lf("init_from_molecule"),
        lf("init_from_rdkit_mol"),
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
