import pytest
import stk

from ._test_case import _TestCase


@pytest.fixture(
    params=(
        stk.N(0),
    )
)
def nitrogen(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.H(1),
    ),
)
def hydrogen1(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.H(2),
    ),
)
def hydrogen2(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.C(3),
    ),
)
def carbon1(request):
    return request.parm.clone()


@pytest.fixture(
    params=(
        stk.C(4),
    ),
)
def carbon2(request):
    return request.param.clone()


@pytest.fixture(
    params=(
        stk.H(5),
    ),
)
def hydrogen3(request):
    return request.param.clone()



@pytest.fixture(
    params=(
        stk.C(6),
    ),
)
def carbon3(request):
    return request.param.clone()


@pytest.fixture
def ring_amine(
    nitrogen,
    hydrogen1,
    hydrogen2,
    carbon1,
    carbon2,
    hydrogen3,
    carbon3,
):
    bonders = ()
    deleters = ()
    return _TestCase(
        functional_group=stk.RingAmine(
            nitrogen=nitrogen,
            hydrogen1=hydrogen1,
            hydrogen2=hydrogen2,
            carbon1=carbon1,
            carbon2=carbon2,
            hydrogen3=hydrogen3,
            carbon3=carbon3,
        ),
        atoms=(
            nitrogen,
            hydrogen1,
            hydrogen2,
            carbon1,
            carbon2,
            hydrogen3,
            carbon3,
        ),
        bonders=bonders,
        deleters=deleters,
    )
