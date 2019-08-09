import pytest
import stk


@pytest.fixture(scope='session')
def amine2():
    return stk.BuildingBlock('NCCCN', ['amine'])


@pytest.fixture
def tmp_amine2():
    return stk.BuildingBlock('NCCCN', ['amine'])


@pytest.fixture(scope='session')
def amine2_conf1():
    amine2 = stk.BuildingBlock('NCCCN', ['amine'])
    amine2.set_position_matrix(amine2.get_position_matrix()*3)
    return amine2


@pytest.fixture(scope='session')
def amine2_alt1():
    return stk.BuildingBlock('NCNCN', ['amine'])


@pytest.fixture(scope='session')
def amine2_alt2():
    return stk.BuildingBlock('NC[Si]CN', ['amine'])


@pytest.fixture(scope='session')
def amine2_alt3():
    return stk.BuildingBlock('NC(Cl)CN', ['amine'])


@pytest.fixture(scope='session')
def amine3():
    return stk.BuildingBlock('NCC(CN)CN', ['amine'])


@pytest.fixture(scope='session')
def amine4():
    return stk.BuildingBlock('NCC(CN)(CN)CN', ['amine'])


@pytest.fixture
def tmp_amine4():
    return stk.BuildingBlock('NCC(CN)(CN)CN', ['amine'])
