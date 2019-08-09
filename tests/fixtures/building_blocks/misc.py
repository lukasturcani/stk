import pytest
import stk
from os.path import join


@pytest.fixture(scope='session')
def boronic_acid2():
    return stk.BuildingBlock('OB(O)c1ccc(B(O)O)nc1', ['boronic_acid'])


@pytest.fixture(scope='session')
def boronic_acid4():
    return stk.BuildingBlock.init_from_file(
        path=join('..', 'data', 'boronic_acid.sdf')
    )


@pytest.fixture(scope='session')
def diol2():
    return stk.BuildingBlock('Oc1cc2cc(O)c(O)nc2cc1O', ['diol'])


@pytest.fixture(scope='session')
def ring_amine():
    return stk.BuildingBlock(
        smiles='Nc1ccc2cc3cc(N)ccc3cc2c1',
        functional_groups=['ring_amine']
    )


@pytest.fixture(scope='session')
def cycle():
    return stk.BuildingBlock('CCCC1CCCCCCCCC1')


@pytest.fixture(scope='session')
def c60():
    return stk.BuildingBlock.init_from_file(
        path=join('..', 'data', 'c60.pdb')
    )


@pytest.fixture(scope='session')
def chained_c60():
    return stk.BuildingBlock.init_from_file(
        path=join('..', 'data', 'chained_c60.mol')
    )
