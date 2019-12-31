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
def oxygen(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.N(2),
    ),
)
def nitrogen(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.H(3),
    )
)
def hydrogen1(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.H(4),
    ),
)
def hydrogen2(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.C(5),
    ),
)
def r(request):
    return request.param.clone()


@pytest.fixture
def amide(carbon, oxygen, nitrogen, hydrogen1, hydrogen2, r):
    bonders = (carbon, )
    deleters = (nitrogen, hydrogen1, hydrogen2)
    return _TestCase(
        functional_group=stk.Amide(
            carbon=carbon,
            oxygen=oxygen,
            nitrogen=nitrogen,
            hydrogen1=hydrogen1,
            hydrogen2=hydrogen2,
            r=r,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(carbon, oxygen, nitrogen, hydrogen1, hydrogen2, r),
        bonders=bonders,
        deleters=deleters,
    )
