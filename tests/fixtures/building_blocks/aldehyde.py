import pytest
import stk


@pytest.fixture(scope='session')
def aldehyde2():
    return stk.BuildingBlock('O=CCC=O', ['aldehyde'])


@pytest.fixture(scope='function')
def tmp_aldehyde2():
    return stk.BuildingBlock('O=CCC=O', ['aldehyde'])


@pytest.fixture('session')
def aldehyde2_alt1():
    return stk.BuildingBlock('O=CCNCC=O', ['aldehyde'])


@pytest.fixture('session')
def aldehyde2_alt2():
    return stk.BuildingBlock('O=CCOCC=O', ['aldehyde'])


@pytest.fixture(scope='session')
def aldehyde3():
    return stk.BuildingBlock('O=CC(C=O)C=O', ['aldehyde'])


@pytest.fixture(scope='function')
def tmp_aldehyde3():
    return stk.BuildingBlock('O=CC(C=O)C=O', ['aldehyde'])


@pytest.fixture(scope='session')
def aldehyde3_alt1():
    return stk.BuildingBlock('O=CN(C=O)C=O', ['aldehyde'])


@pytest.fixture(scope='session')
def aldehyde3_alt2():
    return stk.BuildingBlock('O=C[SiH](C=O)C=O', ['aldehyde'])


@pytest.fixture(scope='session')
def aldehyde3_alt3():
    return stk.BuildingBlock('O=CC(Cl)C(C=O)C=O', ['aldehyde'])


@pytest.fixture(scope='session')
def aldehyde4():
    return stk.BuildingBlock('O=CC(C=O)(C=O)C=O', ['aldehyde'])


@pytest.fixture(scope='function')
def tmp_aldehyde4():
    return stk.BuildingBlock('O=CC(C=O)(C=O)C=O', ['aldehyde'])


@pytest.fixture(scope='session')
def aldehyde4_alt1():
    return stk.BuildingBlock('O=CC(OC=O)(C=O)C=O', ['aldehyde'])


@pytest.fixture(scope='function')
def tmp_aldehyde4_alt2():
    return stk.BuildingBlock(
        smiles='O=Cc1c(C=O)cc(C=O)c(C=O)c1',
        functional_groups=['aldehyde']
    )


@pytest.fixture(scope='session')
def aldehyde5():
    return stk.BuildingBlock(
        smiles='O=C[C-]1C(C=O)=C(C=O)C(C=O)=C1C=O',
        functional_groups=['aldehyde']
    )


@pytest.fixture(scope='function')
def tmp_aldehyde5():
    return stk.BuildingBlock(
        smiles='O=C[C-]1C(C=O)=C(C=O)C(C=O)=C1C=O',
        functional_groups=['aldehyde']
    )


@pytest.fixture(scope='session')
def aldehyde6():
    return stk.BuildingBlock(
        smiles='O=Cc1c(C=O)c(C=O)c(C=O)c(C=O)c1C=O',
        functional_groups=['aldehyde']
    )


@pytest.fixture(scope='function')
def tmp_aldehyde6():
    return stk.BuildingBlock(
        smiles='O=Cc1c(C=O)c(C=O)c(C=O)c(C=O)c1C=O',
        functional_groups=['aldehyde']
    )
