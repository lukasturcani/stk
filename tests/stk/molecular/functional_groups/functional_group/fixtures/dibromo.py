import pytest
import stk

from ._test_case import _TestCase


@pytest.fixture(
    params=(
        stk.Br(0),
    ),
)
def bromine1(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.C(9),
    ),
)
def atom1(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.Br(1),
    ),
)
def bromine2(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.C(21),
    ),
)
def atom2(request):
    return request.param.clone()


@pytest.fixture
def dibromo(bromine1, atom1, bromine2, atom2):
    bonders = (atom1, atom2)
    deleters = (bromine1, bromine2)
    return _TestCase(
        functional_group=stk.Dibromo(
            bromine1=bromine1,
            atom1=atom1,
            bromine2=bromine2,
            atom2=atom2,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(bromine1, atom1, bromine2, atom2),
        bonders=bonders,
        deleters=deleters,
    )
