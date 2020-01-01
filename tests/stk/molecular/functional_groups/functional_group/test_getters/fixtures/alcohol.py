import pytest
import stk

from ._test_case import _TestCase


@pytest.fixture(
    params=(
        stk.O(0),
    ),
)
def oxygen(request):
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
def alcohol(oxygen, hydrogen, atom):
    bonders = (oxygen, )
    deleters = (hydrogen, )
    return _TestCase(
        functional_group=stk.Alcohol(
            oxygen=oxygen,
            hydrogen=hydrogen,
            atom=atom,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(oxygen, hydrogen, atom),
        bonders=bonders,
        deleters=deleters,
    )
