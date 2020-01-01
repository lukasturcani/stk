import pytest
import stk

from ._test_case import _TestCase


@pytest.fixture(
    params=(
        stk.B(0),
    ),
)
def boron(request):
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
        stk.O(3),
    ),
)
def oxygen2(request):
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
def atom(request):
    return request.param.clone()


@pytest.fixture
def boronic_acid(boron, oxygen1, hydrogen1, oxygen2, hydrogen2, atom):
    bonders = (oxygen1, oxygen2)
    deleters = (hydrogen1, hydrogen2)
    return _TestCase(
        functional_group=stk.BoronicAcid(
            boron=boron,
            oxygen1=oxygen1,
            hydrogen1=hydrogen1,
            oxygen2=oxygen2,
            hydrogen2=hydrogen2,
            atom=atom,
        ),
        atoms=(boron, oxygen1, hydrogen1, oxygen2, hydrogen2, atom),
        bonders=bonders,
        deleters=deleters,
    )
