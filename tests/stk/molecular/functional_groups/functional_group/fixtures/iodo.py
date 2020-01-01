import pytest
import stk

from ._test_case import _TestCase


@pytest.fixture(
    params=(
        stk.I(0),
    ),
)
def iodine(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.C(1),
    ),
)
def atom(request):
    return request.param.clone()


@pytest.fixture
def iodo(iodine, atom):
    bonders = (atom, )
    deleters = (iodine, )
    return _TestCase(
        functional_group=stk.Iodo(
            iodine=iodine,
            atom=atom,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(iodine, atom),
        bonders=bonders,
        deleters=deleters,
    )
