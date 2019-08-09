import pytest
import stk


@pytest.fixture(scope='session')
def ab_chain3():
    return stk.polymer.Linear('AB', [0, 0], 3)


@pytest.fixture
def tmp_ab_chain3():
    return stk.polymer.Linear('AB', [0, 0], 3)


@pytest.fixture(scope='session')
def ab_chain6():
    return stk.polymer.Linear('AB', [0, 0], 6)
