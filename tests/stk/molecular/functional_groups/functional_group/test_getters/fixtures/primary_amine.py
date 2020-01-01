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
def hydrogen1(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.H(2),
    )
)
def hydrogen2(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.C(3),
    )
)
def atom(request):
    return request.param.clone()


@pytest.fixture
def primary_amine(nitrogen, hydrogen1, hydrogen2, atom):
    bonders = (nitrogen, )
    deleters = (hydrogen1, hydrogen2)
    return _TestCase(
        functional_group=stk.PrimaryAmine(
            nitrogen=nitrogen,
            hydrogen1=hydrogen1,
            hydrogen2=hydrogen2,
            atom=atom,
        ),
        atoms=(nitrogen, hydrogen1, hydrogen2, atom),
        bonders=bonders,
        deleters=deleters,
    )
