import stk


def test_random_building_block(polymer,
                               amine2,
                               aldehyde2,
                               aldehyde2_alt1):

    mutator = stk.RandomBuildingBlock(
            building_blocks=[aldehyde2_alt1],
            key=lambda mol: mol.func_group_infos[0].name == 'aldehyde'
    )
    mutant = mutator.mutate(polymer)

    assert mutant.__class__ == polymer.__class__
    assert mutant.topology.__class__ == polymer.topology.__class__
    assert any(
        bb.same(aldehyde2_alt1) for bb in mutant.building_blocks
    )
    assert all(
        not bb.same(aldehyde2) for bb in mutant.building_blocks
    )


def test_similar_bb(polymer,
                    amine2,
                    aldehyde2,
                    aldehyde2_alt1,
                    aldehyde2_alt2):

    mutator = stk.SimilarBuildingBlock(
            building_blocks=[aldehyde2_alt1, aldehyde2_alt2],
            duplicate_building_blocks=False,
            key=lambda mol: mol.func_group_infos[0].name == 'aldehyde'
    )
    mutant = mutator.mutate(polymer)

    assert mutant.__class__ == polymer.__class__
    assert mutant.topology.__class__ == polymer.topology.__class__
    assert any(
        bb.same(aldehyde2_alt1) for bb in mutant.building_blocks
    )
    assert all(
        not bb.same(aldehyde2) for bb in mutant.building_blocks
    )

    mutant = mutator.mutate(polymer)

    assert mutant.__class__ == polymer.__class__
    assert mutant.topology.__class__ == polymer.topology.__class__
    assert any(
        bb.same(aldehyde2_alt2) for bb in mutant.building_blocks
    )
    assert all(
        not bb.same(aldehyde2) for bb in mutant.building_blocks
    )


def test_random_topology(cage):
    topologies = [stk.EightPlusTwelve(), stk.FourPlusSix()]
    mutator = stk.RandomTopology(topologies)
    mutant = mutator.mutate(cage)

    assert type(mutant) == type(cage)
    assert mutant.topology.__class__ != cage.topology.__class__
    assert all(
        any(m.same(c) for c in cage.building_blocks)
        for m in mutant.building_blocks
    )
    assert mutant.building_blocks == cage.building_blocks


def test_random_mutation(cage):
    t1 = stk.RandomTopology([stk.EightPlusTwelve()])
    t2 = stk.RandomTopology([stk.FourPlusSix2()])
    mutator1 = stk.RandomMutation(t1, t2, weights=[1, 0])
    mutator2 = stk.RandomMutation(t1, t2, weights=[0, 1])

    mutant1 = mutator1.mutate(cage)
    assert isinstance(mutant1.topology, stk.EightPlusTwelve)
    mutant2 = mutator2.mutate(cage)
    assert isinstance(mutant2.topology, stk.FourPlusSix2)
