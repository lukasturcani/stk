import pytest
import stk
import os
from os.path import join
import logging
import sys
from collections import Counter, defaultdict

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


class TestFitnessCalculator(stk.FitnessCalculator):
    def fitness(self, mol, conformer=-1):
        return 12


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
def fitness_calculator():
    return TestFitnessCalculator(use_cache=False)


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


@pytest.fixture
def tmp_amine2_alt1():
    return stk.StructUnit2.smiles_init(smiles='NCNCN',
                                       functional_groups=['amine'],
                                       name='tmp_amine2_alt1')


@pytest.fixture(scope='session')
def amine2_alt1():
    return stk.BuildingBlock('NCNCN', ['amine'])


@pytest.fixture
def tmp_amine2_alt2():
    return stk.StructUnit2.smiles_init(smiles='NC[Si]CN',
                                       functional_groups=['amine'],
                                       name='tmp_amine2_alt2')


@pytest.fixture(scope='session')
def amine2_alt2():
    return stk.BuildingBlock('NC[Si]CN', ['amine'])


@pytest.fixture(scope='session')
def amine2_alt3():
    return stk.BuildingBlock('NC(Cl)CN', ['amine'])


@pytest.fixture(scope='session')
def aldehyde2():
    return stk.BuildingBlock('O=CCC=O', ['aldehyde'])


@pytest.fixture
def tmp_aldehyde2():
    return stk.BuildingBlock('O=CCC=O', ['aldehyde'])


@pytest.fixture(scope='session')
def aldehyde2_alt1():
    return stk.StructUnit2.smiles_init(smiles='O=CNC=O',
                                       functional_groups=['aldehyde'],
                                       name='aldehyde2_alt1')


@pytest.fixture(scope='session')
def aldehyde2_alt2():
    return stk.StructUnit2.smiles_init(smiles='O=CNNC=O',
                                       functional_groups=['aldehyde'],
                                       name='aldehyde2_alt2')


@pytest.fixture(scope='session')
def boronic_acid2():
    return stk.BuildingBlock('OB(O)c1ccc(B(O)O)nc1', ['boronic_acid'])


@pytest.fixture(scope='session')
def bromine2():
    return stk.StructUnit2.smiles_init(smiles='[Br]CCC[Br]',
                                       functional_groups=['bromine'],
                                       name='bromine2')


@pytest.fixture
def tmp_bromine2():
    return stk.BuildingBlock('[Br]CCC[Br]', ['bromine'])


@pytest.fixture(scope='session')
def bromine2_alt1():
    return stk.StructUnit2.smiles_init(smiles='[Br]CNC[Br]',
                                       functional_groups=['bromine'],
                                       name='bromine2_alt1')


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


@pytest.fixture(scope='session')
def polymer(amine2, aldehyde2):
    return stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2],
        topology_graph=stk.polymer.Linear('AB', [0, 0], 3)
    )


@pytest.fixture(scope='session')
def polymer_alt1(amine2_alt1, aldehyde2_alt1):
    return stk.Polymer([amine2_alt1, aldehyde2_alt1],
                       stk.Linear('AB', [0, 0], 3),
                       'polymer_alt1')


@pytest.fixture(scope='session')
def polymer_alt2(amine2_alt2, aldehyde2_alt2):
    return stk.Polymer([amine2_alt2, aldehyde2_alt2],
                       stk.Linear('AB', [0, 0], 3),
                       'polymer_alt2')


@pytest.fixture
def tmp_polymer(tmp_amine2, tmp_aldehyde2):
    return stk.ConstructedMolecule(
        building_blocks=[tmp_amine2, tmp_aldehyde2],
        topology_graph=stk.polymer.Linear('AB', [0, 0], 3)
    )


@pytest.fixture
def reactor(amine2, aldehyde2):
    building_blocks = [amine2, aldehyde2]
    mol = stk.ConstructedMolecule.__new__(stk.ConstructedMolecule)
    mol.topology_graph = stk.polymer.Linear('AB', [0, 0], 3)
    mol.atoms = []
    mol.bonds = []
    mol.bonds_made = 0
    mol.func_groups = []
    mol.building_block_counter = Counter()
    mol._position_matrix = []
    mol.building_block_vertices = defaultdict(list)
    mol.topology_graph._assign_building_blocks_to_vertices(
        mol=mol,
        building_blocks=building_blocks
    )
    vertex_clones = mol.topology_graph._clone_vertices()
    edge_clones = mol.topology_graph._clone_edges(vertex_clones)
    mol._edge_clones = edge_clones

    mol.topology_graph._prepare(mol)
    mol.topology_graph._place_building_blocks(mol, vertex_clones)
    return stk.molecular.reactor.Reactor(mol)


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


<<<<<<< HEAD
@pytest.fixture(scope='session')
def generate_population(struct_units2, struct_units3):

    def inner(cache=False, offset=False):
        """
        Returns a population of subpopulations and direct members.

        """

        stk.OPTIONS['cache'] = cache

        # Generate a bunch of cages.
        if offset:
            mols = [TestMol([struct_units2[i], struct_units3[i+1]],
                            stk.FourPlusSix(),
                            f'test_mol_{i}')
                    for i in range(10)]

        else:
            mols = [TestMol([struct_units2[i], struct_units3[i]],
                            stk.FourPlusSix(),
                            f'test_mol_{i}')
                    for i in range(10)]

        # Generate a couple of
        # populations to be used as subpopulations.
        sub1 = stk.Population(*mols[:2])
        sub2 = stk.Population(*mols[2:4])
        sub3 = stk.Population(*mols[4:6])
        sub4 = stk.Population(*mols[6:8])

        # Place subpopulations into one another.
        sub1.populations.append(sub3)
        sub2.populations.append(sub4)

        # Initialize final population of subpopulations and cages.
        p = stk.Population(sub1, sub2, *mols[8:])
        p.assign_names_from(0, True)

        stk.OPTIONS['cache'] = False

        return p

    return inner


@pytest.fixture(scope='session')
def pop(generate_population):
    return generate_population()


@pytest.fixture
def tmp_polymer_pop(bromine2, bromine2_alt1):
    p1 = stk.Polymer([bromine2, bromine2_alt1],
                     stk.Linear('AB', [0, 0], 1),
                     'p1')
    p2 = stk.Polymer([bromine2, bromine2_alt1],
                     stk.Linear('ABBA', [0, 0], 1),
                     'p2')
    p3 = stk.Polymer([bromine2, bromine2_alt1],
                     stk.Linear('ABA', [0, 0], 1),
                     'p3')
    p4 = stk.Polymer([bromine2, bromine2_alt1],
                     stk.Linear('AAB', [0, 0], 1),
                     'p4')

    return stk.Population(p1, p2, p3, p4)


@pytest.fixture(scope='session')
def ga_input():
    return stk.GAInput(join('..', 'data', 'inputfile.py'))


@pytest.fixture(scope='session')
def progress():
    pop = stk.Population(*(
        stk.Population(*(
            stk.Molecule(f'{i}') for i in range(20)
        ))
        for i in range(25)
    ))

    for i, mol in enumerate(pop):
        mol.fitness = i
        mol.cavity_size = -i

    return pop


@pytest.fixture(scope='session')
def flat_pop(progress):
    return progress.populations[0]
=======
# @pytest.fixture(scope='session')
# def ga_input():
#     return stk.GAInput(join('..', 'data', 'inputfile.py'))
>>>>>>> fixing_tests
