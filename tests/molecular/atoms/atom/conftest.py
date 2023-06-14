import pytest
import stk
from pytest_lazyfixture import lazy_fixture

from .case_data import CaseData
from .utilities import atomic_numbers


@pytest.fixture(
    params=[
        cls
        for cls in stk.__dict__.values()
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
        atom=stk.Atom(
            id=id,
            atomic_number=atomic_number,
            charge=charge,
        ),
        id=id,
        charge=charge,
        atomic_number=atomic_number,
    )


@pytest.fixture
def case_data_2(cls, id, charge):
    """
    A :class:`.CaseData` instance.

    """

    return CaseData(
        atom=cls(id, charge),
        id=id,
        charge=charge,
        atomic_number=atomic_numbers[cls],
    )


@pytest.fixture(
    params=(
        lazy_fixture("case_data_1"),
        lazy_fixture("case_data_2"),
    )
)
def case_data(request):
    """
    A :class:`.CaseData` instance.

    """

    return request.param


@pytest.fixture
def atom(cls):
    """
    An :class:`.Atom` instance.

    """

    return cls(3, -5).clone()


@pytest.fixture(
    params=[0, 3],
)
def id(request):
    """
    An atom id.

    """

    return request.param


@pytest.fixture(
    params=[0],
)
def charge(request):
    """
    An atomic charge.

    """

    return request.param


@pytest.fixture(params=tuple(range(1, 118)))
def atomic_number(request):
    """
    An atomic number.

    """

    return request.param
