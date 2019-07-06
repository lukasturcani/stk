import pytest
import stk
import os
from os.path import join
import logging
import sys

logging.basicConfig(
    format='\n\n%(levelname)s:%(module)s:%(message)s',
    stream=sys.stdout
)
logging.getLogger('stk').setLevel(logging.DEBUG)


# Run tests in a directory so that that generated files are easy to
# delete.
output_dir = 'tests_output'
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
os.chdir(output_dir)


def pytest_addoption(parser):
    parser.addoption('--macromodel_path', default='')
    parser.addoption('--mopac_path', default='')
    parser.addoption('--xtb_path', default='')


def pytest_generate_tests(metafunc):
    if 'macromodel_path' in metafunc.fixturenames:
        mm_path = metafunc.config.getoption('macromodel_path')
        metafunc.parametrize('macromodel_path', [mm_path])

    if 'mopac_path' in metafunc.fixturenames:
        mopac_path = metafunc.config.getoption('mopac_path')
        metafunc.parametrize('mopac_path', [mopac_path])

    if 'xtb_path' in metafunc.fixturenames:
        xtb_path = metafunc.config.getoption('xtb_path')
        metafunc.parametrize('xtb_path', [xtb_path])


@pytest.fixture(scope='session')
def amine2():
    amine2 = stk.BuildingBlock.init_from_smiles(
        smiles='NCCCN',
        functional_groups=['amine']
    )
    # Make a second conformer with a distinct geometry.
    amine2.set_position_matrix(
        position_matrix=amine2.get_position_matrix()*4,
        conformer_id=None
    )
    return amine2


@pytest.fixture
def tmp_amine2():
    amine2 = stk.BuildingBlock.init_from_smiles(
        smiles='NCCCN',
        functional_groups=['amine']
    )
    # Make a second conformer with a distinct geometry.
    amine2.set_position_from_matrix(
        position_matrix=amine2.get_position_matrix()*4,
        conformer_id=None
    )
    return amine2


@pytest.fixture(scope='session')
def amine2_alt1():
    return stk.BuildingBlock.init_from_smiles(
        smiles='NCNCN',
        functional_groups=['amine']
    )


@pytest.fixture(scope='session')
def amine2_alt2():
    return stk.BuildingBlock.init_from_smiles(
        smiles='NC[Si]CN',
        functional_groups=['amine']
    )


@pytest.fixture(scope='session')
def amine2_alt3():
    return stk.BuildingBlock.init_from_smiles(
        smiles='NC(Cl)CN',
        functional_groups=['amine']
    )


@pytest.fixture(scope='session')
def aldehyde2():
    return stk.BuildingBlock.init_from_smiles(
        smiles='O=CCC=O',
        functional_groups=['aldehyde']
    )


@pytest.fixture(scope='session')
def boronic_acid2():
    return stk.BuildingBlock.init_from_smiles(
        smiles='OB(O)c1ccc(B(O)O)nc1',
        functional_groups=['boronic_acid']
    )


@pytest.fixture
def tmp_bromine2():
    return stk.BuildingBlock.init_from_smiles(
        smiles='[Br]CCC[Br]',
        functional_groups=['bromine']
    )


@pytest.fixture
def tmp_bromine2_alt1():
    return stk.BuildingBlock.init_from_smiles(
        smiles='[Br]CNC[Br]',
        functional_groups=['bromine']
    )


@pytest.fixture(scope='session')
def diol2():
    return stk.BuildingBlock.init_from_smiles(
        smiles='Oc1cc2cc(O)c(O)nc2cc1O',
        functional_groups=['diol']
    )


@pytest.fixture(scope='session')
def amine3():
    return stk.BuildingBlock.init_from_smiles(
        smiles='NCC(CN)CN',
        functional_groups=['amine']
    )


@pytest.fixture(scope='session')
def ring_amine():
    smiles = 'Nc1ccc2cc3cc(N)ccc3cc2c1'
    return stk.BuildingBlock.init_from_smiles(
        smiles=smiles,
        functional_groups=['ring_amine']
    )


@pytest.fixture(scope='session')
def aldehyde3():
    aldehyde3 = stk.BuildingBlock.init_from_smiles(
        smiles='O=CC(C=O)C=O',
        functional_groups=['aldehyde']
    )
    # Make a second conformer with a distinct geometry.
    aldehyde3.set_position_from_matrix(
        position_matrix=aldehyde3.get_position_matrix()*4,
        conformer_id=None
    )
    return aldehyde3


@pytest.fixture
def tmp_aldehyde3():
    aldehyde3 = stk.BuildingBlock.init_from_smiles(
        smiles='O=CC(C=O)C=O',
        functional_groups=['aldehyde']
    )
    # Make a second conformer with a distinct geometry.
    aldehyde3.set_position_from_matrix(
        position_matrix=aldehyde3.get_position_matrix()*4,
        conformer_id=None
    )
    return aldehyde3


@pytest.fixture(scope='session')
def aldehyde3_alt1():
    return stk.BuildingBlock.init_from_smiles(
        smiles='O=CN(C=O)C=O',
        functional_groups=['aldehyde']
    )


@pytest.fixture(scope='session')
def aldehyde3_alt2():
    return stk.BuildingBlock.init_from_smiles(
        smiles='O=C[Si](C=O)C=O',
        functional_groups=['aldehyde']
    )


@pytest.fixture(scope='session')
def aldehyde3_alt3():
    return stk.BuildingBlock.init_from_smiles(
        smiles='O=CC(Cl)C(C=O)C=O',
        functional_groups=['aldehyde']
    )


@pytest.fixture(scope='session')
def boronic_acid4():
    return stk.BuildingBlock(
        mol=join('..', 'data', 'boronic_acid.sdf')
    )


@pytest.fixture(scope='session')
def amine4():
    return stk.BuildingBlock.init_from_smiles(
        smiles='NCC(CN)(CN)CN',
        functional_groups=['amine']
    )


@pytest.fixture(scope='session')
def aldehyde4():
    return stk.BuildingBlock.init_from_smiles(
        smiles='O=CC(C=O)(C=O)C=O',
        functional_groups=['aldehyde']
    )


@pytest.fixture(scope='session')
def aldehyde4_alt1():
    return stk.BuildingBlock.init_from_smiles(
        smiles='O=CC(OC=O)(C=O)C=O',
        functional_groups=['aldehyde']
    )


@pytest.fixture(scope='session')
def difluorene2():
    return stk.BuildingBlock.init_from_smiles(
        smiles='Fc1c(F)cc(F)c(F)c1',
        functional_groups=['difluorene']
    )


@pytest.fixture(scope='session')
def difluorene_dibromine():
    return stk.BuildingBlock.init_from_smiles(
        smiles='Fc1c(F)cc(Br)c(Br)c1',
        functional_groups=['difluorene', 'dibromine']
    )


@pytest.fixture(scope='session')
def aldehyde5():
    return stk.BuildingBlock.init_from_smiles(
        smiles='O=C[C-]1C(C=O)=C(C=O)C(C=O)=C1C=O',
        functional_groups=['aldehyde']
    )


@pytest.fixture(scope='session')
def aldehyde6():
    return stk.BuildingBlock.init_from_smiles(
        smiles='O=CC(C=O)(C=O)C(C=O)(C=O)C=O',
        functional_groups=['aldehyde']
    )


@pytest.fixture(scope='session')
def cycle():
    return stk.BuildingBlock.init_from_smiles(
        smiles='CCCC1CCCCCCCCC1'
    )


@pytest.fixture(scope='session')
def polymer(amine2, aldehyde2):
    return stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2],
        topology=stk.polymer.Linear('AB', [0, 0], 3)
    )


@pytest.fixture
def tmp_polymer(amine2, aldehyde2):
    return stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2],
        topology=stk.polymer.Linear('AB', [0, 0], 3)
    )


@pytest.fixture(scope='session')
def cc3():
    # This has an unoptimized conformer in conformer 0 and an
    # optimized one in conformer 1.
    m = stk.Molecule.load(join('..', 'data', 'cc3.json'))
    m.name = 'cc3'
    return m


@pytest.fixture
def tmp_cc3():
    # This has an unoptimized conformer in conformer 0 and an
    # optimized one in conformer 1.
    m = stk.Molecule.load(join('..', 'data', 'cc3.json'))
    m.name = 'tmp_cc3'
    return m


@pytest.fixture(scope='session')
def c60():
    return stk.BuildingBlock(join('..', 'data', 'c60.pdb'))


@pytest.fixture(scope='session')
def chained_c60():
    return stk.BuildingBlock(join('..', 'data', 'chained_c60.mol'))


@pytest.fixture(scope='session')
def fg():
    return stk.FunctionalGroup(
        id=0,
        atom_ids=[10, 3, 1, 4, 43, 5, 32, 55],
        bonder_ids=[3, 32, 10],
        deleter_ids=[1, 55, 5],
        info=stk.functional_groups[0]
    )


@pytest.fixture
def tmp_fg():
    return stk.FunctionalGroup(
        id=0,
        atom_ids=[10, 3, 1, 4, 43, 5, 32, 55],
        bonder_ids=[3, 32, 10],
        deleter_ids=[1, 55, 5],
        info=stk.functional_groups[0]
    )


@pytest.fixture(scope='session')
def mae_path():
    return join('..', 'data', 'molecule.mae')


@pytest.fixture(scope='session')
def population():
    return stk.Population(*(
        stk.BuildingBlock.init_from_smiles('C'*i)
        for i in range(1, 12)
    ))


@pytest.fixture(scope='session')
def ga_input():
    return stk.GAInput(join('..', 'data', 'inputfile.py'))
