import pytest
from pytest_lazyfixture import lazy_fixture
import stk

from .utilities import atomic_numbers, atomic_masses
from .case_data import CaseData


@pytest.fixture(
    params=[
        cls for cls in stk.__dict__.values()
        if isinstance(cls, type)
        and issubclass(cls, stk.Atom)
        and cls is not stk.Atom
    ],
)
def cls(request):
    """
    Return an :class:`.Atom` instance.

    Parameters
    ----------
    id : :class:`int`
        The id of the returned atom.

    charge : :class:`float`
        The charge of the returned atom.

    Returns
    -------
    :class:`.Atom`
        An atom.

    """

    return request.param


@pytest.fixture
def case_data_1(atomic_number, id, charge):
    """
    A :class:`.CaseData` instance.

    """

    return CaseData(
        atom=stk.Atom(atomic_number, id, charge),
        id=id,
        charge=charge,
        atomic_number=atomic_number,
        mass=atomic_masses[atomic_number],
    )


def case_data_2(cls, id, charge):
    return CaseData(
        atom=cls(id, charge),
        id=id,
        charge=charge,
        atomic_number=atomic_numbers[cls],
        mass=atomic_masses[atomic_numbers[cls]],
    )


@pytest.fixture(
    params=(
        lazy_fixture('case_data_1'),
        lazy_fixture('case_data_2'),
    )
)
def case_data(request):
    return request.param


@pytest.fixture
def atom(get_atom):
    """
    An :class:`.Atom` instance.

    """

    return get_atom(3, -5).clone()


@pytest.fixture(
    params=[0, 10],
)
def id(request):
    """
    An atom id.

    """

    return request.param


@pytest.fixture(
    params=[0, 10, -1],
)
def charge(request):
    """
    An atomic charge.

    """

    return request.param
