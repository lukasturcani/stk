import pytest
import stk

from ._test_case import _TestCase


@pytest.fixture(
    params=(
        stk.C(0),
    ),
)
def atom(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.F(1),
    ),
)
def fluorine(request):
    return request.param.clone()


@pytest.fixture
def fluoro(fluorine, atom):
    bonders = (atom, )
    deleters = (fluorine, )
    return _TestCase(
        functional_group=stk.Fluoro(
            fluorine=fluorine,
            atom=atom,
            bonders=bonders,
            deleters=deleters,
        ),
        atoms=(fluorine, atom),
        bonders=bonders,
        deleters=deleters,
    )
