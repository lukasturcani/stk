import pytest
from pytest_lazy_fixtures import lf

# Fixtures need to be visible for lazy_fixture() calls.
from .fixtures import *


@pytest.fixture(
    params=(
        lf("primary_amino"),
        lf("secondary_amino"),
        lf("single_atom"),
        lf("aldehyde"),
        lf("carboxylic_acid"),
        lf("amide"),
        lf("thioacid"),
        lf("alcohol"),
        lf("alkene"),
        lf("alkyne"),
        lf("boronic_acid"),
        lf("bromo"),
        lf("dibromo"),
        lf("difluoro"),
        lf("diol"),
        lf("fluoro"),
        lf("iodo"),
        lf("thiol"),
    ),
)
def generic_case_data(request):
    """
    A :class:`.GenericCaseData` instance.

    """

    return request.param


@pytest.fixture(
    params=(
        lf("generic_case_data"),
        lf("ring_amine"),
    ),
)
def case_data(request):
    """
    A :class:`.CaseData` instance.

    """

    return request.param


@pytest.fixture(
    params=(lambda n: range(n),),
)
def get_atom_ids(request):
    """
    Yield `n` atom ids.

    Parameters
    ----------
    n : :class:`int`
        The number of atom ids to yield.

    Yields
    ------
    :class:`int`
        An atom id.

    """

    return request.param


@pytest.fixture
def functional_group(case_data):
    """
    A :class:`.FunctionalGroup` instance.

    """

    return case_data.functional_group


@pytest.fixture
def generic_functional_group(generic_case_data):
    """
    A :class:`.GenericFunctionalGroup` instance.

    """

    return generic_case_data.functional_group
