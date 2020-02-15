import pytest

from .fixtures import *


@pytest.fixture(
    params=(
        pytest.lazy_fixture('primary_amino'),
        pytest.lazy_fixture('secondary_amino'),
        pytest.lazy_fixture('aldehyde'),
        pytest.lazy_fixture('carboxylic_acid'),
        pytest.lazy_fixture('amide'),
        pytest.lazy_fixture('thioacid'),
        pytest.lazy_fixture('alcohol'),
        pytest.lazy_fixture('alkene'),
        pytest.lazy_fixture('alkyne'),
        pytest.lazy_fixture('boronic_acid'),
        pytest.lazy_fixture('bromo'),
        pytest.lazy_fixture('dibromo'),
        pytest.lazy_fixture('difluoro'),
        pytest.lazy_fixture('diol'),
        pytest.lazy_fixture('fluoro'),
        pytest.lazy_fixture('iodo'),
        pytest.lazy_fixture('thiol'),
    ),
)
def generic_test_case(request):
    """
    A :class:`._GenericTestCase` instance.

    """

    return request.param


@pytest.fixture(
    params=(
        pytest.lazy_fixture('generic_test_case'),
        pytest.lazy_fixture('ring_amine'),
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
