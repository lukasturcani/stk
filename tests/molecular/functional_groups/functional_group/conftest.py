import pytest
from pytest_lazyfixture import lazy_fixture

# Fixtures need to be visible for lazy_fixture() calls.
from .fixtures import * # noqa


@pytest.fixture(
    params=(
        lazy_fixture('primary_amino'),
        lazy_fixture('secondary_amino'),
        lazy_fixture('aldehyde'),
        lazy_fixture('carboxylic_acid'),
        lazy_fixture('amide'),
        lazy_fixture('thioacid'),
        lazy_fixture('alcohol'),
        lazy_fixture('alkene'),
        lazy_fixture('alkyne'),
        lazy_fixture('boronic_acid'),
        lazy_fixture('bromo'),
        lazy_fixture('dibromo'),
        lazy_fixture('difluoro'),
        lazy_fixture('diol'),
        lazy_fixture('fluoro'),
        lazy_fixture('iodo'),
        lazy_fixture('thiol'),
    ),
)
def generic_test_case(request):
    """
    A :class:`._GenericTestCase` instance.

    """

    return request.param


@pytest.fixture(
    params=(
        lazy_fixture('generic_test_case'),
        lazy_fixture('ring_amine'),
    ),
)
def test_case(request):
    """
    A :class:`._TestCase` instance.

    """

    return request.param


@pytest.fixture(
    params=(
        lambda n: range(n),
    ),
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
def functional_group(test_case):
    """
    A :class:`.FunctionalGroup` instance.

    """

    return test_case.functional_group


@pytest.fixture
def generic_functional_group(generic_test_case):
    """
    A :class:`.GenericFunctionalGroup` instance.

    """

    return generic_test_case.functional_group
