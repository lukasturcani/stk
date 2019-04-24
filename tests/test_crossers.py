import stk


def test_random_crossover(polymer, polymer_alt1):
    jumble = stk.Jumble(num_offspring_building_blocks=2)
    recombination = stk.GeneticRecombination(
                        key=lambda mol: mol.func_group_infos[0].name)
    random_crossover = stk.RandomCrossover(jumble,
                                           recombination,
                                           weights=[0, 1])
    cohort1 = list(random_crossover.crossover(polymer, polymer_alt1))
    assert len(cohort1) == 4

    random_crossover = stk.RandomCrossover(jumble,
                                           recombination,
                                           weights=[1, 0])
    cohort1 = list(random_crossover.crossover(polymer, polymer_alt1))
    assert len(cohort1) == 6


def test_jumble(polymer, polymer_alt1, polymer_alt2):
    jumble = stk.Jumble(num_offspring_building_blocks=2)
    cohort1 = list(jumble.crossover(polymer, polymer_alt1))

    assert len(cohort1) == 6
    for mol1 in cohort1:
        assert len(mol1.building_blocks) == 2
        for mol2 in cohort1:
            if mol1 is not mol2:
                assert not mol1.same(mol2)

    jumble = stk.Jumble(num_offspring_building_blocks=2,
                        duplicate_building_blocks=True)
    cohort2 = list(jumble.crossover(polymer, polymer_alt1))

    assert len(cohort2) == 10
    for mol1 in cohort2:
        assert len(mol1.building_blocks) == 2
        for mol2 in cohort2:
            if mol1 is not mol2:
                assert not mol1.same(mol2)

    jumble = stk.Jumble(num_offspring_building_blocks=2)
    cohort3 = list(
        jumble.crossover(polymer, polymer_alt1, polymer_alt2)
    )

    assert len(cohort3) == 15
    for mol1 in cohort3:
        assert len(mol1.building_blocks) == 2
        for mol2 in cohort3:
            if mol1 is not mol2:
                assert not mol1.same(mol2)


def test_genetic_recombination(polymer, polymer_alt1, polymer_alt2):
    recombination = stk.GeneticRecombination(
                        key=lambda mol: mol.func_group_infos[0].name)

    cohort1 = list(recombination.crossover(polymer, polymer_alt1))

    assert len(cohort1) == 4
    for mol1 in cohort1:
        assert len(mol1.building_blocks) == 2
        for mol2 in cohort1:
            if mol1 is not mol2:
                assert not mol1.same(mol2)

    cohort2 = list(
        recombination.crossover(polymer, polymer_alt1, polymer_alt2)
    )

    assert len(cohort2) == 9
    for mol1 in cohort2:
        assert len(mol1.building_blocks) == 2
        for mol2 in cohort2:
            if mol1 is not mol2:
                assert not mol1.same(mol2)
