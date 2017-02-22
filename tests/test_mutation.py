from ..ga import Mutation
from ..molecular import Molecule, EightPlusTwelve, FourPlusSix
from os.path import join

path = join('data', 'mutation', 'm.json')
with open(path, 'r') as f:
    mol = Molecule.load(f.read())

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
