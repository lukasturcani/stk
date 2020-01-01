import pytest
import stk

from ._test_case import _TestCase


@pytest.fixture(
    params=(
        stk.S(0),
    ),
)
def sulfur(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.H(1),
    ),
)
def hydrogen(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.C(2),
    ),
)
def atom(request):
    return request.param.clone()


@pytest.fixture
def thiol(sulfur, hydrogen, atom):
    bonders = ()
    deleters = ()
    return _TestCase(
        functional_group=stk.Thiol(
            sulfur=sulfur,
            hydrogen=hydrogen,
            atom=atom,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(sulfur, hydrogen, atom),
        bonders=bonders,
        deleters=deleters,
    )
