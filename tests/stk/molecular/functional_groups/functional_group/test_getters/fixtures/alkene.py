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
        stk.H(1),
    ),
)
def atom1(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.C(2),
    ),
)
def atom2(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.C(3),
    ),
)
def carbon2(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.H(4),
    ),
)
def atom3(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.C(5),
    ),
)
def atom4(request):
    return request.param.clone()


@pytest.fixture
def alkene(carbon1, atom1, atom2, carbon2, atom3, atom4):
    bonders = (carbon2, )
    deleters = (carbon1, atom1, atom2)
    return _TestCase(
        functional_group=stk.Alkene(
            carbon1=carbon1,
            atom1=atom1,
            atom2=atom2,
            carbon2=carbon2,
            atom3=atom3,
            atom4=atom4,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(carbon1, atom1, atom2, carbon2, atom3, atom4),
        bonders=bonders,
        deleters=deleters,
    )
