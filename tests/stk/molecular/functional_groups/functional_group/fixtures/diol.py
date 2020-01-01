import pytest
import stk

from ._test_case import _TestCase


@pytest.fixture(
    params=(
        stk.C(0),
    ),
)
def atom1(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.O(1),
    ),
)
def oxygen1(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.H(2),
    ),
)
def hydrogen1(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.C(3),
    ),
)
def atom2(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.O(4),
    ),
)
def oxygen2(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.H(5),
    ),
)
def hydrogen2(request):
    return request.param.clone()


@pytest.fixture
def diol(atom1, oxygen1, hydrogen1, atom2, oxygen2, hydrogen2):
    bonders = (atom1, atom2)
    deleters = (oxygen1, hydrogen1, oxygen2, hydrogen2)
    return _TestCase(
        functional_group=stk.Diol(
            atom1=atom1,
            oxygen1=oxygen1,
            hydrogen1=hydrogen1,
            atom2=atom2,
            oxygen2=oxygen2,
            hydrogen2=hydrogen2,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(atom1, oxygen1, hydrogen1, atom2, oxygen2, hydrogen2),
        bonders=bonders,
        deleters=deleters,
    )
