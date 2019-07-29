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


@pytest.fixture('session')
def hydrogen():
    return stk.H(1, -3, attr1='12', attr2=32, _attr3=12.2)


@pytest.fixture('session')
def carbon():
    return stk.C(333, attr2=None, attr12='abc', _pattr=22.22)


@pytest.fixture('session')
def lithium():
    return stk.Li(121, 2, alpha=232, beta='bvc', _gamma=4523)


@pytest.fixture('session')
def chlorine():
    return stk.Cl(786, a=9, b='bbb', _c='yfg')


@pytest.fixture('session')
def bond(hydrogen, carbon):
    return stk.Bond(
        atom1=hydrogen,
        atom2=carbon,
        order=2,
        attr1=1,
        attr2='2',
        _attr3=12.2
    )


@pytest.fixture('session')
def periodic_bond(lithium, chlorine):
    return stk.PeriodicBond(
        atom1=lithium,
        atom2=chlorine,
        order=21,
        direction=[1, 0, -1],
        attr10=16,
        attr20='26',
        _attr30=126.2
    )


@pytest.fixture(scope='session')
def amine2():
    return stk.BuildingBlock('NCCCN', ['amine'])


@pytest.fixture
def tmp_amine2():
    return stk.BuildingBlock('NCCCN', ['amine'])


@pytest.fixture
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
def aldehyde2():
    return stk.BuildingBlock('O=CCC=O', ['aldehyde'])


@pytest.fixture(scope='session')
def boronic_acid2():
    return stk.BuildingBlock('OB(O)c1ccc(B(O)O)nc1', ['boronic_acid'])


@pytest.fixture
def tmp_bromine2():
    return stk.BuildingBlock('[Br]CCC[Br]', ['bromine'])


@pytest.fixture
def tmp_bromine2_alt1():
    return stk.BuildingBlock('[Br]CNC[Br]', ['bromine'])


@pytest.fixture(scope='session')
def diol2():
    return stk.BuildingBlock('Oc1cc2cc(O)c(O)nc2cc1O', ['diol'])


@pytest.fixture(scope='session')
def amine3():
    return stk.BuildingBlock('NCC(CN)CN', ['amine'])


@pytest.fixture(scope='session')
def ring_amine():
    return stk.BuildingBlock.init_from_smiles(
        smiles='Nc1ccc2cc3cc(N)ccc3cc2c1',
        functional_groups=['ring_amine']
    )


@pytest.fixture(scope='session')
def aldehyde3():
    return stk.BuildingBlock('O=CC(C=O)C=O', ['aldehyde'])


@pytest.fixture
def tmp_aldehyde3():
    return stk.BuildingBlock('O=CC(C=O)C=O', ['aldehyde'])


@pytest.fixture(scope='session')
def aldehyde3_alt1():
    return stk.BuildingBlock('O=CN(C=O)C=O', ['aldehyde'])


@pytest.fixture(scope='session')
def aldehyde3_alt2():
    return stk.BuildingBlock('O=C[Si](C=O)C=O', ['aldehyde'])


@pytest.fixture(scope='session')
def aldehyde3_alt3():
    return stk.BuildingBlock('O=CC(Cl)C(C=O)C=O', ['aldehyde'])


@pytest.fixture(scope='session')
def boronic_acid4():
    return stk.BuildingBlock.init_from_file(
        path=join('..', 'data', 'boronic_acid.sdf')
    )


@pytest.fixture(scope='session')
def amine4():
    return stk.BuildingBlock('NCC(CN)(CN)CN', ['amine'])


@pytest.fixture
def tmp_amine4():
    return stk.BuildingBlock('NCC(CN)(CN)CN', ['amine'])


@pytest.fixture(scope='session')
def aldehyde4():
    return stk.BuildingBlock('O=CC(C=O)(C=O)C=O', ['aldehyde'])


@pytest.fixture(scope='session')
def aldehyde4_alt1():
    return stk.BuildingBlock('O=CC(OC=O)(C=O)C=O', ['aldehyde'])


@pytest.fixture(scope='session')
def difluorene2():
    return stk.BuildingBlock('Fc1c(F)cc(F)c(F)c1', ['difluorene'])


@pytest.fixture(scope='session')
def difluorene_dibromine():
    return stk.BuildingBlock(
        smiles='Fc1c(F)cc(Br)c(Br)c1',
        functional_groups=['difluorene', 'dibromine']
    )


@pytest.fixture(scope='session')
def aldehyde5():
    return stk.BuildingBlock(
        smiles='O=C[C-]1C(C=O)=C(C=O)C(C=O)=C1C=O',
        functional_groups=['aldehyde']
    )


@pytest.fixture(scope='session')
def aldehyde6():
    return stk.BuildingBlock(
        smiles='O=CC(C=O)(C=O)C(C=O)(C=O)C=O',
        functional_groups=['aldehyde']
    )


@pytest.fixture(scope='session')
def cycle():
    return stk.BuildingBlock('CCCC1CCCCCCCCC1')


# @pytest.fixture(scope='session')
# def polymer(amine2, aldehyde2):
#     return stk.ConstructedMolecule(
#         building_blocks=[amine2, aldehyde2],
#         topology=stk.polymer.Linear('AB', [0, 0], 3)
#     )
#
#
# @pytest.fixture
# def tmp_polymer(amine2, aldehyde2):
#     return stk.ConstructedMolecule(
#         building_blocks=[amine2, aldehyde2],
#         topology=stk.polymer.Linear('AB', [0, 0], 3)
#     )


# @pytest.fixture(scope='session')
# def cc3():
#     m = stk.Molecule.load(join('..', 'data', 'cc3.json'))
#     m.name = 'cc3'
#     return m


# @pytest.fixture
# def tmp_cc3():
#     m = stk.Molecule.load(join('..', 'data', 'cc3.json'))
#     m.name = 'tmp_cc3'
#     return m


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


@pytest.fixture(scope='session')
def mae_path():
    return join('..', 'data', 'molecule.mae')


@pytest.fixture(scope='session')
def bb_dir():
    return join('..', 'data', 'building_block_init')


# @pytest.fixture(scope='session')
# def population():
#     return stk.Population(*(
#         stk.BuildingBlock.init_from_smiles('C'*i)
#         for i in range(1, 12)
#     ))


# @pytest.fixture(scope='session')
# def ga_input():
#     return stk.GAInput(join('..', 'data', 'inputfile.py'))
