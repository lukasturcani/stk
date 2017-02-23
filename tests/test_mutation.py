from ..ga import Mutation, Population
from ..molecular import EightPlusTwelve, FourPlusSix, StructUnit
from os.path import join
from glob import glob
path = join('data', 'mutation', 'mutants.json')
pop = Population.load(path)
mol = pop[0]



def test_cage_random_bb():

    for bb_file in glob(r'data/mutation/cage/bb/*.mol'):
        StructUnit(bb_file, 'aldehyde')
        

    mutant = Mutation.cage_random_bb(None, mol,
                join('data', 'mutation', 'cage', 'bb'),
                'aldehyde')

    _, bb1 = min(zip(mol.bb_counter.values(),
                     mol.bb_counter.keys()))
    _, lk1 = max(zip(mol.bb_counter.values(),
                     mol.bb_counter.keys()))

    assert mutant.topology.__class__ == mol.topology.__class__
    assert False

def test_cage_random_lk():
    mutant = Mutation.cage_random_lk(None, mol,
                join('data', 'mutation', 'cage', 'lk'),
                     'amine')

    _, bb1 = min(zip(mol.bb_counter.values(),
                     mol.bb_counter.keys()))
    _, lk1 = max(zip(mol.bb_counter.values(),
                     mol.bb_counter.keys()))

    assert mutant.topology.__class__ == mol.topology.__class__
    assert False

def test_cage_similar_bb():
    ...

def test_cage_similar_lk():
    ...

def test_random_topology():
    mutant = Mutation.random_topology(None, mol, [EightPlusTwelve(),
                                                  FourPlusSix()])
    assert type(mutant) == type(mol)
    assert mutant.topology.__class__ != mol.topology.__class__
    assert mutant.building_blocks == mol.building_blocks
