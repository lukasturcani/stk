import pytest
import stk
import os
from os.path import join
from collections import Counter
import rdkit.Chem.AllChem as rdkit
import logging
import sys

logging.basicConfig(format='\n\n%(levelname)s:%(module)s:%(message)s',
                    stream=sys.stdout)
logging.getLogger('stk').setLevel(logging.DEBUG)


# Run tests in a directory so that that generated files are easy to
# delete.
output_dir = 'tests_output'
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
os.chdir(output_dir)

stk.OPTIONS['cache'] = False


class TestMol(stk.Cage):
    def __init__(self,
                 building_blocks,
                 topology,
                 name="",
                 note="",
                 bb_conformers=None):
        if bb_conformers is None:
            bb_conformers = [-1 for _ in range(len(building_blocks))]

        self.building_blocks = building_blocks
        self.bb_counter = Counter(building_blocks)
        self.bonds_made = len(building_blocks) - 1
        self.mol = rdkit.Mol()
        for bb in building_blocks:
            self.mol = rdkit.CombineMols(self.mol, bb.mol)
        self.topology = topology
        stk.Molecule.__init__(self, name, note)
        self.func_groups = (
            stk.FunctionalGroup(id_=0,
                                atom_ids=[0, 1, 2, 3],
                                bonder_ids=[0, 1],
                                deleter_ids=[2],
                                info='amine'),
            stk.FunctionalGroup(id_=1,
                                atom_ids=[10, 20, 30],
                                bonder_ids=[10],
                                deleter_ids=[20, 30],
                                info='aldehyde')
        )

    def windows(self, *_, **__):
        return [4, 3, 2, 1]

    def cavity_size(self, *_, **__):
        return 3.48

    def bb_distortion(self, *_, **__):
        return 5

    def dihedral_strain(self, *_, **__):
        return 4


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


@pytest.fixture(scope='session')
def fitness_calculator():
    return TestFitnessCalculator(use_cache=False)


@pytest.fixture(scope='session')
def amine2():
    amine2 = stk.StructUnit2.smiles_init(smiles='NCCCN',
                                         functional_groups=['amine'],
                                         name='amine2')
    # Make a second conformer with a distinct geometry.
    amine2.mol.AddConformer(amine2.mol.GetConformer(), True)
    amine2.set_position_from_matrix(
        pos_mat=amine2.mol.GetConformer().GetPositions().T*4,
        conformer=1
    )
    return amine2


@pytest.fixture
def tmp_amine2():
    amine2 = stk.StructUnit2.smiles_init(smiles='NCCCN',
                                         functional_groups=['amine'],
                                         name='tmp_amine2')
    # Make a second conformer with a distinct geometry.
    amine2.mol.AddConformer(amine2.mol.GetConformer(), True)
    amine2.set_position_from_matrix(
        pos_mat=amine2.mol.GetConformer().GetPositions().T*4,
        conformer=1
    )
    return amine2


@pytest.fixture
def tmp_amine2_alt1():
    return stk.StructUnit2.smiles_init(smiles='NCNCN',
                                       functional_groups=['amine'],
                                       name='tmp_amine2_alt1')


@pytest.fixture(scope='session')
def amine2_alt1():
    return stk.StructUnit2.smiles_init(smiles='NCNCN',
                                       functional_groups=['amine'],
                                       name='amine2_alt1')


@pytest.fixture
def tmp_amine2_alt2():
    return stk.StructUnit2.smiles_init(smiles='NC[Si]CN',
                                       functional_groups=['amine'],
                                       name='tmp_amine2_alt2')


@pytest.fixture(scope='session')
def amine2_alt2():
    return stk.StructUnit2.smiles_init(smiles='NC[Si]CN',
                                       functional_groups=['amine'],
                                       name='amine2_alt2')


@pytest.fixture(scope='session')
def amine2_alt3():
    return stk.StructUnit2.smiles_init(smiles='NC(Cl)CN',
                                       functional_groups=['amine'],
                                       name='amine2_alt3')


@pytest.fixture(scope='session')
def aldehyde2():
    return stk.StructUnit2.smiles_init(smiles='O=CCC=O',
                                       functional_groups=['aldehyde'],
                                       name='aldehyde2')


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
    return stk.StructUnit2.smiles_init(
                                smiles='OB(O)c1ccc(B(O)O)nc1',
                                functional_groups=['boronic_acid'],
                                name='boronic_acid2')


@pytest.fixture(scope='session')
def bromine2():
    return stk.StructUnit2.smiles_init(smiles='[Br]CCC[Br]',
                                       functional_groups=['bromine'],
                                       name='bromine2')


@pytest.fixture
def tmp_bromine2():
    return stk.StructUnit2.smiles_init(smiles='[Br]CCC[Br]',
                                       functional_groups=['bromine'],
                                       name='tmp_bromine2')


@pytest.fixture(scope='session')
def bromine2_alt1():
    return stk.StructUnit2.smiles_init(smiles='[Br]CNC[Br]',
                                       functional_groups=['bromine'],
                                       name='bromine2_alt1')


@pytest.fixture
def tmp_bromine2_alt1():
    return stk.StructUnit2.smiles_init(smiles='[Br]CNC[Br]',
                                       functional_groups=['bromine'],
                                       name='tmp_bromine2_alt1')


@pytest.fixture(scope='session')
def diol2():
    return stk.StructUnit2.smiles_init(smiles='Oc1cc2cc(O)c(O)nc2cc1O',
                                       functional_groups=['diol'],
                                       name='diol2')


@pytest.fixture(scope='session')
def amine3():
    return stk.StructUnit3.smiles_init(smiles='NCC(CN)CN',
                                       functional_groups=['amine'],
                                       name='amine3')


@pytest.fixture(scope='session')
def ring_amine():
    smiles = 'Nc1ccc2cc3cc(N)ccc3cc2c1'
    return stk.StructUnit2.smiles_init(
                                smiles=smiles,
                                functional_groups=['ring_amine'],
                                name='ring_amine')


@pytest.fixture(scope='session')
def aldehyde3():
    aldehyde3 = stk.StructUnit3.smiles_init(
                                smiles='O=CC(C=O)C=O',
                                functional_groups=['aldehyde'],
                                name='aldehyde3')
    # Make a second conformer with a distinct geometry.
    aldehyde3.mol.AddConformer(aldehyde3.mol.GetConformer(), True)
    aldehyde3.set_position_from_matrix(
        pos_mat=aldehyde3.mol.GetConformer().GetPositions().T*4,
        conformer=1
    )
    return aldehyde3


@pytest.fixture
def tmp_aldehyde3():
    aldehyde3 = stk.StructUnit3.smiles_init(
                                smiles='O=CC(C=O)C=O',
                                functional_groups=['aldehyde'],
                                name='tmp_aldehyde3')
    # Make a second conformer with a distinct geometry.
    aldehyde3.mol.AddConformer(aldehyde3.mol.GetConformer(), True)
    aldehyde3.set_position_from_matrix(
        pos_mat=aldehyde3.mol.GetConformer().GetPositions().T*4,
        conformer=1
    )
    return aldehyde3


@pytest.fixture(scope='session')
def aldehyde3_alt1():
    return stk.StructUnit3.smiles_init(smiles='O=CN(C=O)C=O',
                                       functional_groups=['aldehyde'],
                                       name='aldehyde3_alt1')


@pytest.fixture(scope='session')
def aldehyde3_alt2():
    return stk.StructUnit3.smiles_init(smiles='O=C[Si](C=O)C=O',
                                       functional_groups=['aldehyde'],
                                       name='aldehyde3_alt2')


@pytest.fixture(scope='session')
def aldehyde3_alt3():
    return stk.StructUnit3.smiles_init(smiles='O=CC(Cl)C(C=O)C=O',
                                       functional_groups=['aldehyde'],
                                       name='aldehyde3_alt3')


@pytest.fixture(scope='session')
def boronic_acid4():
    return stk.StructUnit3(join('..', 'data', 'boronic_acid.sdf'),
                           name='boronic_acid4')


@pytest.fixture(scope='session')
def amine4():
    return stk.StructUnit3.smiles_init(smiles='NCC(CN)(CN)CN',
                                       functional_groups=['amine'],
                                       name='amine4')


@pytest.fixture(scope='session')
def aldehyde4():
    return stk.StructUnit3.smiles_init(smiles='O=CC(C=O)(C=O)C=O',
                                       functional_groups=['aldehyde'],
                                       name='aldehyde4')


@pytest.fixture(scope='session')
def aldehyde4_alt1():
    return stk.StructUnit3.smiles_init(smiles='O=CC(OC=O)(C=O)C=O',
                                       functional_groups=['aldehyde'],
                                       name='aldehyde4_alt1')


@pytest.fixture(scope='session')
def difluorene2():
    smiles = 'Fc1c(F)cc(F)c(F)c1'
    return stk.StructUnit2.smiles_init(
                                    smiles=smiles,
                                    functional_groups=['difluorene'],
                                    name='difluorene2')


@pytest.fixture(scope='session')
def difluorene_dibromine():
    smiles = 'Fc1c(F)cc(Br)c(Br)c1'
    return stk.StructUnit2.smiles_init(
                    smiles=smiles,
                    functional_groups=['difluorene', 'dibromine'],
                    name='difluorene_dibromine')


@pytest.fixture(scope='session')
def cage(amine2, aldehyde3):
    return stk.Cage([amine2, aldehyde3], stk.FourPlusSix(), 'cage')


@pytest.fixture
def tmp_cage(amine2, aldehyde3):
    return stk.Cage([amine2, aldehyde3], stk.FourPlusSix(), 'tmp_cage')


@pytest.fixture(scope='session')
def aldehyde6():
    return stk.StructUnit3.smiles_init(
                                smiles='O=CC(C=O)(C=O)C(C=O)(C=O)C=O',
                                functional_groups=['aldehyde'],
                                name='aldehyde6')


@pytest.fixture(scope='session')
def cycle_su():
    cycle = rdkit.MolFromSmiles('CCCC1CCCCCCCCC1')
    cycle = rdkit.AddHs(cycle)
    rdkit.EmbedMolecule(cycle)

    return stk.MacrocycleStructUnit(cycle, [])


@pytest.fixture(scope='session')
def cycle(amine2, aldehyde2):
    return stk.Macrocycle([amine2, aldehyde2],
                          stk.Cyclic('AB', [0, 0], 3))


@pytest.fixture
def tmp_cycle(amine2, aldehyde2):
    return stk.Macrocycle([amine2, aldehyde2],
                          stk.Cyclic('AB', [0, 0], 3))


@pytest.fixture(scope='session')
def polymer(amine2, aldehyde2):
    return stk.Polymer([amine2, aldehyde2],
                       stk.Linear('AB', [0, 0], 3),
                       'polymer')


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
def tmp_polymer(amine2, aldehyde2):
    return stk.Polymer([amine2, aldehyde2],
                       stk.Linear('AB', [0, 0], 3),
                       'tmp_polymer')


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
    return stk.StructUnit(join('..', 'data', 'c60.pdb'), name='c60')


@pytest.fixture(scope='session')
def chained_c60():
    path = join('..', 'data', 'chained_c60.mol')
    return stk.StructUnit(path, name='c60')


@pytest.fixture(scope='session')
def fg():
    return stk.FunctionalGroup(id_=0,
                               atom_ids=[10, 3, 1, 4, 43, 5, 32, 55],
                               bonder_ids=[3, 32, 10],
                               deleter_ids=[1, 55, 5],
                               info=stk.functional_groups[0])


@pytest.fixture
def tmp_fg():
    return stk.FunctionalGroup(id_=0,
                               atom_ids=[10, 3, 1, 4, 43, 5, 32, 55],
                               bonder_ids=[3, 32, 10],
                               deleter_ids=[1, 55, 5],
                               info=stk.functional_groups[0])


@pytest.fixture(scope='session')
def test_mol1():
    bb1 = stk.StructUnit2.smiles_init('N')
    bb2 = stk.StructUnit3.smiles_init('NN')
    # Make sure calling build does nothing.
    top = stk.FourPlusSix()
    top.build = lambda x: ...
    test_mol = TestMol([bb1, bb2], top, 'test_mol1')
    test_mol.mol = rdkit.Mol(bb1.mol)
    test_mol.bonds_made = 2
    return test_mol


@pytest.fixture(scope='session')
def mae_path():
    return join('..', 'data', 'molecule.mae')


@pytest.fixture(scope='session')
def struct_units2():
    return [
        stk.StructUnit2.smiles_init('C'*i, name=f'{i}')
        for i in range(1, 12)
    ]


@pytest.fixture(scope='session')
def struct_units3():
    return [
        stk.StructUnit3.smiles_init('C'*i, name=f'{i}')
        for i in range(1, 12)
    ]


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
