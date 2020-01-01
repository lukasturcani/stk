import pytest
import stk

from ._test_case import _TestCase


@pytest.fixture(
    params=(
        stk.Br(0),
    )
)
def bromine(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.C(1),
    ),
)
def atom(request):
    return request.param.clone()


@pytest.fixture
def bromo(bromine, atom):
    bonders = (atom, )
    deleters = (bromine, )
    return _TestCase(
        functional_group=stk.Bromo(
            bromine=bromine,
            atom=atom,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(bromine, atom),
        bonders=bonders,
        deleters=deleters,
    )
