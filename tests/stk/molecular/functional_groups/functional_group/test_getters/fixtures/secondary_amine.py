import pytest
import stk

from ._test_case import _TestCase


@pytest.fixture(
    params=(
        stk.N(0),
    ),
)
def nitrogen(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.H(1),
    )
)
def hydrogen(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.C(2),
    )
)
def r1(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.C(3),
    )
)
def r2(request):
    return request.param.clone()


@pytest.fixture
def secondary_amine(nitrogen, hydrogen, r1, r2):
    bonders = (nitrogen, )
    deleters = (hydrogen, )
    return _TestCase(
        functional_group=stk.SecondaryAmine(
            nitrogen=nitrogen,
            hydrogen=hydrogen,
            r1=r1,
            r2=r2,
        ),
        atoms=(nitrogen, hydrogen, r1, r2),
        bonders=bonders,
        deleters=deleters,
    )
