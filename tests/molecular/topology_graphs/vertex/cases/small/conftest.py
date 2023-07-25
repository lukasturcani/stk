import numpy as np
import pytest
import stk
from pytest_lazyfixture import lazy_fixture

from ...case_data import CaseData


@pytest.fixture(
    params=(
        lazy_fixture("substituent"),
        lazy_fixture("core"),
    ),
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param


@pytest.fixture
def core(
    id: int,
    position: tuple[float, float, float] | np.ndarray,
) -> CaseData:
    return CaseData(
        vertex=stk.small.CoreVertex(
            id=id,
            position=position,
        ),
        id=id,
        position=position,
        cell=np.array([0, 0, 0]),
    )


@pytest.fixture
def substituent(
    id: int,
    position: tuple[float, float, float] | np.ndarray,
) -> CaseData:
    return CaseData(
        vertex=stk.small.SubstituentVertex(
            id=id,
            position=position,
        ),
        id=id,
        position=position,
        cell=np.array([0, 0, 0]),
    )


@pytest.fixture(params=(0,))
def id(request: pytest.FixtureRequest) -> int:
    return request.param
