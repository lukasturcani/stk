from ..ga import Mutation, Population
from ..molecular import (EightPlusTwelve,
                         FourPlusSix,
                         Molecule,
                         StructUnit3)
from os.path import join
from glob import iglob

path = join('data', 'mutation', 'mutants.json')
pop = Population.load(path, Molecule.from_dict)
mol = pop[0]


def test_random_bb():
    mols = [StructUnit3(path, 'aldehyde') for path in
            iglob(join('data', 'mutation', 'cage', 'bb', '*.mol'))]
    mutant = Mutation(None, None).random_bb(
                                mol,
                                mols,
                                lambda x: x.__class__ is StructUnit3)

    assert mutant.topology.__class__ == mol.topology.__class__


def test_similar_bb():
    db = join('data', 'mutation', 'cage', 'bb')
    mols = [StructUnit3(path, 'aldehyde') for
            path in iglob(join(db, '*.mol'))]

    mutant = Mutation(None, None).similar_bb(
                                mol,
                                mols,
                                lambda x: x.__class__ is StructUnit3)

    rdkit_mols = (m.mol for m in mols)
    most_sim = mol.building_blocks[0].similar_molecules(rdkit_mols)[1][1]
    mol_map = {m.mol: m for m in mols}
    assert mol_map[most_sim] is mutant.building_blocks[1]
    assert mutant.topology.__class__ == mol.topology.__class__


def test_random_topology():
    mutant = Mutation.random_topology(None, mol, [EightPlusTwelve(),
                                                  FourPlusSix()])
    assert type(mutant) == type(mol)
    assert mutant.topology.__class__ != mol.topology.__class__
    assert mutant.building_blocks == mol.building_blocks
