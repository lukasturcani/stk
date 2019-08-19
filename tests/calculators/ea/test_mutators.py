import stk


def test_random_building_block(
    polymer,
    amine2,
    aldehyde2,
    aldehyde2_alt1
):
    mutator = stk.RandomBuildingBlock(
        building_blocks=[aldehyde2_alt1],
        key=lambda mol: mol.func_groups[0].fg_type.name == 'aldehyde'
    )
    mutant = mutator.mutate(polymer)
    expected = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2_alt1],
        topology_graph=stk.polymer.Linear('AB', 3)
    )
    assert mutant.get_identity_key() == expected.get_identity_key()


def test_similar_building_block(
    polymer,
    amine2,
    aldehyde2,
    aldehyde2_alt1,
    aldehyde2_alt2
):

    mutator = stk.SimilarBuildingBlock(
        building_blocks=[aldehyde2_alt1, aldehyde2_alt2],
        duplicate_building_blocks=False,
        key=lambda mol: mol.func_groups[0].fg_type.name == 'aldehyde'
    )
    mutant = mutator.mutate(polymer)
    expected = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2_alt1],
        topology_graph=stk.polymer.Linear('AB', 3)
    )
    assert mutant.get_identity_key() == expected.get_identity_key()

    mutant = mutator.mutate(polymer)
    expected = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2_alt2],
        topology_graph=stk.polymer.Linear('AB', 3)
    )
    assert mutant.get_identity_key() == expected.get_identity_key()


def test_random_topology_graph(polymer):
    chain1 = stk.polymer.Linear('AB', 4)
    chain2 = stk.polymer.Linear('AB', 5)
    expected1 = stk.ConstructedMolecule(
        building_blocks=list(polymer.building_block_vertices.keys()),
        topology_graph=chain1
    )
    expected2 = stk.ConstructedMolecule(
        building_blocks=list(polymer.building_block_vertices.keys()),
        topology_graph=chain2
    )
    mutator = stk.RandomTopologyGraph([chain1, chain2])
    mutant = mutator.mutate(polymer)

    expected = {
        expected1.get_identity_key(),
        expected2.get_identity_key()
    }
    assert mutant.get_identity_key() in expected


def test_random_mutation(polymer):
    chain1 = stk.polymer.Linear('AB', 4)
    chain2 = stk.polymer.Linear('AB', 5)
    t1 = stk.RandomTopologyGraph([chain1])
    t2 = stk.RandomTopologyGraph([chain2])
    mutator1 = stk.RandomMutation(t1, t2, weights=[1, 0])
    mutator2 = stk.RandomMutation(t1, t2, weights=[0, 1])

    mutant1 = mutator1.mutate(polymer)
    expected1 = stk.ConstructedMolecule(
        building_blocks=list(polymer.building_block_vertices.keys()),
        topology_graph=chain1
    )
    assert mutant1.get_identity_key() == expected1.get_identity_key()

    expected2 = stk.ConstructedMolecule(
        building_blocks=list(polymer.building_block_vertices.keys()),
        topology_graph=chain2
    )
    mutant2 = mutator2.mutate(polymer)
    assert mutant2.get_identity_key() == expected2.get_identity_key()
