import pytest
import stk


@pytest.fixture
def atom(atom_cls):
    return atom_cls(3, -5).clone()


@pytest.fixture(
    params=[
        cls for cls in stk.__dict__.values()
        if isinstance(cls, type)
        and issubclass(cls, stk.Atom)
        and cls is not stk.Atom
    ],
)
def atom_cls(request):
    return request.param


@pytest.fixture(
    params=[0, 10],
)
def id(request):
    return request.param


@pytest.fixture(
    params=[0, 10, -1],
)
def charge(request):
    return request.param
