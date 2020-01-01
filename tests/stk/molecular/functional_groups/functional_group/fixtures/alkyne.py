import stk
import pytest

from ._test_case import _TestCase


@pytest.fixture(
    params=(
        stk.C(0),
    ),
)
def carbon1(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.C(1),
    ),
)
def atom1(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.C(2),
    ),
)
def carbon2(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.H(3),
    ),
)
def atom2(request):
    return request.param.clone()


@pytest.fixture
def alkyne(carbon1, atom1, carbon2, atom2):
    bonders = (carbon1, )
    deleters = (carbon2, atom2)
    return _TestCase(
        functional_group=stk.Alkyne(
            carbon1=carbon1,
            atom1=atom1,
            carbon2=carbon2,
            atom2=atom2,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(carbon1, atom1, carbon2, atom2),
        bonders=bonders,
        deleters=deleters,
    )
