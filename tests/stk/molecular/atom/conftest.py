import pytest
import stk


@pytest.fixture(
    params=[
        cls for cls in stk.__dict__.values()
        if isinstance(cls, type)
        and issubclass(cls, stk.Atom)
        and cls is not stk.Atom
    ],
)
def get_atom(request):
    """
    A function which returns an :class:`.Atom` instance.

    The function takes two parameters, the id and charge the created
    atom should have.

    """

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
