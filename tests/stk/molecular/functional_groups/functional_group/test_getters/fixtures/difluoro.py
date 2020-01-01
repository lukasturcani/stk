import pytest
import stk

from ._test_case import _TestCase


@pytest.fixture(
    params=(
        stk.F(1),
    ),
)
def fluorine1(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.C(0),
    ),
)
def atom1(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.F(2),
    ),
)
def fluorine2(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.C(3),
    ),
)
def atom2(request):
    return request.param.clone()


@pytest.fixture
def difluoro(fluorine1, atom1, fluorine2, atom2):
    bonders = (atom1, atom2)
    deleters = (fluorine1, fluorine2)
    return _TestCase(
        functional_group=stk.Difluoro(
            fluorine1=fluorine1,
            atom1=atom1,
            fluorine2=fluorine2,
            atom2=atom2,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(fluorine1, atom1, fluorine2, atom2),
        bonders=bonders,
        deleters=deleters,
    )
