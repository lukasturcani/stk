import stk
import pytest

from ._test_case import _TestCase


@pytest.fixture(
    params=(
        stk.C(0),
    ),
)
def carbon(request):
    return request.clone()


@pytest.fixture(
    params=(
        stk.O(1),
    ),
)
def oxygen(request):
    return request.clone()


@pytest.fixture(
    params=(
        stk.H(2),
    ),
)
def hydrogen(request):
    return request.clone()


@pytest.fixture(
    params=(
        stk.C(3),
    ),
)
def r(request):
    return request.clone()


@pytest.fixture
def aldehyde(carbon, oxygen, hydrogen, r):
    bonders = (carbon, )
    deleters = (oxygen, )
    return _TestCase(
        functional_group=stk.Aldehyde(
            carbon=carbon,
            oxygen=oxygen,
            hydrogen=hydrogen,
            r=r,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(carbon, oxygen, hydrogen, r),
        bonders=bonders,
        deleters=deleters,
    )
