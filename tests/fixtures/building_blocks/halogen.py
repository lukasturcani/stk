import pytest
import stk


@pytest.fixture(scope='function')
def tmp_bromine2():
    return stk.BuildingBlock('[Br]CCC[Br]', ['bromine'])


@pytest.fixture(scope='function')
def tmp_bromine2_alt1():
    return stk.BuildingBlock('[Br]CNC[Br]', ['bromine'])


@pytest.fixture(scope='session')
def difluorene2():
    return stk.BuildingBlock('Fc1c(F)cc(F)c(F)c1', ['difluorene'])


@pytest.fixture(scope='session')
def difluorene_dibromine():
    return stk.BuildingBlock(
        smiles='Fc1c(F)cc(Br)c(Br)c1',
        functional_groups=['difluorene', 'dibromine']
    )
