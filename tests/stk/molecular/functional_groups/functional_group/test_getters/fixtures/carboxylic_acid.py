import stk
import pytest

from ._test_case import _TestCase


@pytest.fixture(
    params=(
        stk.C(0),
    ),
)
def carbon(request):
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
        stk.O(2),
    ),
)
def oxygen2(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.H(3),
    ),
)
def hydrogen(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.C(4),
    ),
)
def r(request):
    return request.param.clone()


@pytest.fixture
def carboxylic_acid(carbon, oxygen1, oxygen2, hydrogen, r):
    bonders = (carbon, )
    deleters = (oxygen2, hydrogen)
    return _TestCase(
        functional_group=stk.CarboxylicAcid(
            carbon=carbon,
            oxygen1=oxygen1,
            oxygen2=oxygen2,
            hydrogen=hydrogen,
            r=r,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(carbon, oxygen1, oxygen2, hydrogen),
        bonders=bonders,
        deleters=deleters,
    )
