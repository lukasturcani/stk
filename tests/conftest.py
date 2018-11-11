import pytest
import stk
from os.path import join
from collections import Counter
import rdkit.Chem.AllChem as rdkit

stk.OPTIONS['cache'] = False


class TestMol(stk.MacroMolecule):
    def __init__(self,
                 building_blocks,
                 topology,
                 name="",
                 note="",
                 bb_conformers=None):
        if bb_conformers is None:
            bb_conformers = [-1 for _ in range(len(building_blocks))]

        self.fitness = None
        self.unscaled_fitness = {}
        self.progress_params = {}
        self.building_blocks = building_blocks
        self.bb_counter = Counter(building_blocks)
        self.bonds_made = len(building_blocks) - 1
        self.mol = rdkit.Mol()
        for bb in building_blocks:
            self.mol = rdkit.CombineMols(self.mol, bb.mol)
        self.topology = topology
        stk.Molecule.__init__(self, name, note)


@pytest.fixture(scope='session')
def amine2():
    return stk.StructUnit2.smiles_init('NCCCN', 'amine')


@pytest.fixture
def tmp_amine2():
    return stk.StructUnit2.smiles_init('NCCCN', 'amine')


@pytest.fixture(scope='session')
def amine2_alt1():
    return stk.StructUnit2.smiles_init('NCNCN', 'amine')


@pytest.fixture(scope='session')
def amine2_alt2():
    return stk.StructUnit2.smiles_init('NC[Si]CN', 'amine')


@pytest.fixture(scope='session')
def aldehyde2():
    return stk.StructUnit2.smiles_init('O=CCC=O', 'aldehyde')


@pytest.fixture(scope='session')
def boronic_acid2():
    return stk.StructUnit2.smiles_init('OB(O)c1ccc(B(O)O)nc1',
                                       'boronic_acid')


@pytest.fixture(scope='session')
def diol2():
    return stk.StructUnit2.smiles_init('Oc1cc2cc(O)c(O)nc2cc1O',
                                       'diol')


@pytest.fixture(scope='session')
def amine3():
    return stk.StructUnit3.smiles_init('NCC(CN)CN', 'amine')


@pytest.fixture(scope='session')
def aldehyde3():
    return stk.StructUnit3.smiles_init('O=CC(C=O)C=O', 'aldehyde')


@pytest.fixture
def tmp_aldehyde3():
    return stk.StructUnit3.smiles_init('O=CC(C=O)C=O', 'aldehyde')


@pytest.fixture(scope='session')
def aldehyde3_alt1():
    return stk.StructUnit3.smiles_init('O=CN(C=O)C=O', 'aldehyde')


@pytest.fixture(scope='session')
def aldehyde3_alt2():
    return stk.StructUnit3.smiles_init('O=C[Si](C=O)C=O', 'aldehyde')


@pytest.fixture(scope='session')
def boronic_acid4():
    return stk.StructUnit3(join('data', 'cof', 'boronic_acid.sdf'))


@pytest.fixture(scope='session')
def amine4():
    return stk.StructUnit3.smiles_init('NCC(CN)(CN)CN', 'amine')


@pytest.fixture(scope='session')
def aldehyde4():
    return stk.StructUnit3.smiles_init('O=CC(C=O)(C=O)C=O', 'aldehyde')


@pytest.fixture(scope='session')
def aldehyde6():
    return stk.StructUnit3.smiles_init('O=CC(C=O)(C=O)C(C=O)(C=O)C=O',
                                       'aldehyde')


@pytest.fixture(scope='session')
def polymer(amine2, aldehyde2):
    return stk.Polymer([amine2, aldehyde2],
                       stk.Linear('AB', [0, 0], 1))


@pytest.fixture(scope='session')
def cc3():
    return stk.Molecule.load(join('data', 'cage', 'cc3.json'))


@pytest.fixture(scope='session')
def generate_population():

    def inner(cache=False, offset=False):
        """
        Returns a population of subpopulations and direct members.

        """

        stk.OPTIONS['cache'] = cache

        bbs = [stk.StructUnit.smiles_init('C'*i)
               for i in range(1, 12)]

        # Generate a bunch of cages.
        if offset:
            mols = [TestMol([bbs[i], bbs[i+1]], stk.FourPlusSix())
                    for i in range(0, 10)]

        if not offset:
            mols = [TestMol([bbs[x], bbs[x]], stk.FourPlusSix())
                    for x in range(0, 10)]

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

        stk.OPTIONS['cache'] = False

        return p

    return inner


@pytest.fixture(scope='session')
def pop(pop_generator):
    return generate_population()


@pytest.fixture(scope='session')
def cage():
    return stk.Molecule.load(join('data', 'fitness', 'cage.json'))


@pytest.fixture(scope='session')
def target():
    return join('data', 'fitness', 'target.pdb')


@pytest.fixture(scope='session')
def ga_input():
    return stk.GAInput(join('data', 'gainput', 'test.py'))
