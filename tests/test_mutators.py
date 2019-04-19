import stk


def test_random_bb(test_mol1, struct_units3):
    mutant = stk.Mutation(None, None).random_bb(
                            test_mol1,
                            struct_units3,
                            lambda x: x.__class__ is stk.StructUnit3)

    assert mutant.__class__ == test_mol1.__class__
    assert mutant.topology.__class__ == test_mol1.topology.__class__
    mol_struct_unit2 = next(bb for bb in test_mol1.building_blocks if
                            isinstance(bb, stk.StructUnit2))
    mol_struct_unit3 = next(bb for bb in test_mol1.building_blocks if
                            isinstance(bb, stk.StructUnit3))

    mutant_struct_unit2 = next(bb for bb in mutant.building_blocks if
                               isinstance(bb, stk.StructUnit2))
    mutant_struct_unit3 = next(bb for bb in mutant.building_blocks if
                               isinstance(bb, stk.StructUnit3))

    assert mol_struct_unit2 == mutant_struct_unit2
    assert mol_struct_unit3 != mutant_struct_unit3


def test_similar_bb(test_mol1, struct_units3):
    mutant = stk.Mutation(None, None).similar_bb(
                            test_mol1,
                            struct_units3,
                            lambda x: x.__class__ is stk.StructUnit3)

    rdkit_mols = (m.mol for m in struct_units3)
    bb3 = next(bb for bb in test_mol1.building_blocks if
               isinstance(bb, stk.StructUnit3))
    mutant_bb3 = next(bb for bb in mutant.building_blocks if
                      isinstance(bb, stk.StructUnit3))
    most_sim = bb3.similar_molecules(rdkit_mols)[0][1]
    mol_map = {m.mol: m for m in struct_units3}
    assert mol_map[most_sim] is mutant_bb3
    assert mutant.topology.__class__ == test_mol1.topology.__class__


def test_random_topology(test_mol1):
    topologies = [stk.EightPlusTwelve(), stk.FourPlusSix()]
    mutant = stk.Mutation.random_topology(None, test_mol1, topologies)

    assert type(mutant) == type(test_mol1)
    assert mutant.topology.__class__ != test_mol1.topology.__class__
    assert mutant.building_blocks == test_mol1.building_blocks
