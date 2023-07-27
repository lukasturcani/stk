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
    ),
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param


@pytest.fixture
def center(
    id: int,
    position: tuple[float, float, float] | np.ndarray,
    flip: bool,
) -> CaseData:
    return CaseData(
        vertex=stk.polymer.LinearVertex(id, position, flip),
        id=id,
        position=position,
        cell=np.array([0, 0, 0]),
    )


@pytest.fixture
def head(
    id: int,
    position: tuple[float, float, float] | np.ndarray,
    flip: bool,
) -> CaseData:
    return CaseData(
        vertex=stk.polymer.HeadVertex(id, position, flip),
        id=id,
        position=position,
        cell=np.array([0, 0, 0]),
    )


@pytest.fixture
def tail(
    id: int,
    position: tuple[float, float, float] | np.ndarray,
    flip: bool,
) -> CaseData:
    return CaseData(
        vertex=stk.polymer.TailVertex(id, position, flip),
        id=id,
        position=position,
        cell=np.array([0, 0, 0]),
    )


@pytest.fixture
def unaligning(
    id: int,
    position: tuple[float, float, float] | np.ndarray,
    flip: bool,
) -> CaseData:
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
def id(request: pytest.FixtureRequest) -> int:
    return request.param


@pytest.fixture(params=(True, False))
def flip(request: pytest.FixtureRequest) -> bool:
    return request.param
