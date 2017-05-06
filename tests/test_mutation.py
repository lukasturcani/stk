from ..ga import Mutation, Population
from ..molecular import (EightPlusTwelve, FourPlusSix, StructUnit,
StructUnit2, StructUnit3)
from os.path import join
from glob import glob

path = join('data', 'mutation', 'mutants.json')
pop = Population.load(path)
mol = pop[0]


def test_cage_random_bb():
    mutant = Mutation(None, None).cage_random_bb(mol,
                join('data', 'mutation', 'cage', 'bb'),
                'aldehyde')

    assert mutant.topology.__class__ == mol.topology.__class__


def test_cage_random_lk():
    mutant = Mutation(None, None).cage_random_lk(mol,
                join('data', 'mutation', 'cage', 'lk'),
                     'amine')

    assert mutant.topology.__class__ == mol.topology.__class__


def test_cage_similar_bb():
    db = join('data', 'mutation', 'cage', 'bb')
    mutant = Mutation(None, None).cage_similar_bb(mol, db, 'aldehyde')


    _, bb1 = min(zip(mol.bb_counter.values(),
                     mol.bb_counter.keys()))
    _, lk1 = max(zip(mol.bb_counter.values(),
                     mol.bb_counter.keys()))

    _, mbb1 = min(zip(mutant.bb_counter.values(),
                     mutant.bb_counter.keys()))
    _, mlk1 = max(zip(mutant.bb_counter.values(),
                     mutant.bb_counter.keys()))

    most_sim = bb1.similar_molecules(db)[0][1]
    assert StructUnit3(most_sim, 'aldehyde') is mbb1
    assert mutant.topology.__class__ == mol.topology.__class__


def test_cage_similar_lk():
    db = join('data', 'mutation', 'cage', 'lk')
    mutant = Mutation(None, None).cage_similar_lk(mol, db, 'amine')


    _, bb1 = min(zip(mol.bb_counter.values(),
                     mol.bb_counter.keys()))
    _, lk1 = max(zip(mol.bb_counter.values(),
                     mol.bb_counter.keys()))

    _, mbb1 = min(zip(mutant.bb_counter.values(),
                     mutant.bb_counter.keys()))
    _, mlk1 = max(zip(mutant.bb_counter.values(),
                     mutant.bb_counter.keys()))

    most_sim = bb1.similar_molecules(db)[0][1]
    assert StructUnit2(most_sim, 'amine') is mlk1
    assert mutant.topology.__class__ == mol.topology.__class__


def test_random_topology():
    mutant = Mutation.random_topology(None, mol, [EightPlusTwelve(),
                                                  FourPlusSix()])
    assert type(mutant) == type(mol)
    assert mutant.topology.__class__ != mol.topology.__class__
    assert mutant.building_blocks == mol.building_blocks
