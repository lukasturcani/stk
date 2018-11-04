import pytest
import stk
from os.path import join

stk.CACHE_SETTINGS['ON'] = False


@pytest.fixture(scope='session')
def mol():
    return stk.StructUnit.smiles_init('NC1CC(Br)C(Br)CC1N')


@pytest.fixture
def tmp_mol():
    return stk.StructUnit.smiles_init('NC1CC(Br)C(Br)CC1N')


@pytest.fixture(scope='session')
def mol2():
    stk.StructUnit.smiles_init('NC1CC(N)CC(N)C1', 'amine')


@pytest.fixture(scope='session')
def mol3():
    return stk.StructUnit2.smiles_init('NCCCN', 'amine')


@pytest.fixture(scope='session')
def mol4():
    return stk.StructUnit3.smiles_init('NCC(N)CN', 'amine')


@pytest.fixture(scope='session')
def cof_bb1():
    return stk.StructUnit2.smiles_init('Nc1ccc(N)cc1', 'amine')


@pytest.fixture(scope='session')
def cc3():
    return stk.Molecule.load(join('data', 'molecule', 'cc3.json'))


@pytest.fixture(scope='session')
def pop():
    ...


@pytest.fixture(scope='session')
def struct_unit2_mols():
    ...


@pytest.fixture(scope='session')
def struct_unit3_mols():
    ...
