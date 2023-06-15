import numpy as np
import pytest
import stk
from pytest_lazyfixture import lazy_fixture

from ...case_data import CaseData


@pytest.fixture(
    params=(
        lazy_fixture("center"),
        lazy_fixture("head"),
        lazy_fixture("tail"),
        lazy_fixture("unaligning"),
    ),
)
def case_data(request):
    return request.param


@pytest.fixture
def center(id, position, flip):
    return CaseData(
        vertex=stk.polymer.LinearVertex(id, position, flip),
        id=id,
        position=position,
        cell=np.array([0, 0, 0]),
    )


@pytest.fixture
def head(id, position, flip):
    return CaseData(
        vertex=stk.polymer.HeadVertex(id, position, flip),
        id=id,
        position=position,
        cell=np.array([0, 0, 0]),
    )


@pytest.fixture
def tail(id, position, flip):
    return CaseData(
        vertex=stk.polymer.TailVertex(id, position, flip),
        id=id,
        position=position,
        cell=np.array([0, 0, 0]),
    )


@pytest.fixture
def unaligning(id, position, flip):
    return CaseData(
        vertex=stk.polymer.UnaligningVertex(
            id=id,
            position=position,
            flip=flip,
        ),
        id=id,
        position=position,
        cell=np.array([0, 0, 0]),
    )


@pytest.fixture(params=(0,))
def id(request):
    return request.param


@pytest.fixture(params=(True, False))
def flip(request):
    return request.param
