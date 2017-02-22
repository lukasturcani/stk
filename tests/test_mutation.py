from ..ga import Mutation, Population
from ..molecular import EightPlusTwelve, FourPlusSix
from os.path import join

path = join('data', 'mutation', 'molecules.json')
pop = Population.load(path)
mol = pop[0]

def test_cage_random_bb():
    ...

def test_cage_random_lk():
    ...

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
