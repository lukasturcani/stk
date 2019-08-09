import stk
from collections import Counter


def test_random_crossover(polymer, polymer_alt1):
    jumble = stk.Jumble(num_offspring_building_blocks=2)
    recombination = stk.GeneticRecombination(
        key=lambda mol: mol.func_groups[0].fg_type.name
    )
    random_crossover = stk.RandomCrossover(
        jumble,
        recombination,
        weights=[0, 1]
    )
    cohort1 = list(random_crossover.cross(polymer, polymer_alt1))
    assert len(cohort1) == 6

    random_crossover = stk.RandomCrossover(
        jumble,
        recombination,
        weights=[1, 0]
    )
    cohort1 = list(random_crossover.cross(polymer, polymer_alt1))
    assert len(cohort1) == 10


def _test_cohort(cohort, parents, duplicates):
    topologies = set()
    parents = {parent._key for parent in parents}
    bb_counts = Counter()
    for mol1 in cohort:
        assert mol1._key not in parents
        bb_counts.update([len(mol1.building_block_vertices)])
        topologies.add(repr(mol1.topology_graph))
        for mol2 in cohort:
            if mol1 is not mol2:
                assert not mol1.is_identical(mol2)
    assert len(topologies) == 2

    if duplicates:
        assert len(bb_counts) == 2
        assert bb_counts[1] == 4*len(parents)
        assert bb_counts[2] == len(cohort) - len(parents)*4
    else:
        assert len(bb_counts) == 1
        assert bb_counts[2] == len(cohort)


def test_jumble(polymer, polymer_alt1, polymer_alt2):
    jumble = stk.Jumble(num_offspring_building_blocks=2)
    parents = [polymer, polymer_alt1]

    cohort1 = list(jumble.cross(*parents))
    assert len(cohort1) == 10
    _test_cohort(cohort1, parents, False)

    jumble = stk.Jumble(
        num_offspring_building_blocks=2,
        duplicate_building_blocks=True
    )
    cohort2 = list(jumble.cross(*parents))
    assert len(cohort2) == 18
    _test_cohort(cohort2, parents, True)

    jumble = stk.Jumble(num_offspring_building_blocks=2)
    parents = [polymer, polymer_alt1, polymer_alt2]
    cohort3 = list(jumble.cross(*parents))
    assert len(cohort3) == 27
    _test_cohort(cohort3, parents, False)


def test_genetic_recombination(polymer, polymer_alt1, polymer_alt2):
    recombination = stk.GeneticRecombination(
        key=lambda mol: mol.func_groups[0].fg_type.name
    )
    parents = [polymer, polymer_alt1]
    cohort1 = list(recombination.cross(*parents))
    assert len(cohort1) == 6
    _test_cohort(cohort1, parents, False)

    parents = [polymer, polymer_alt1, polymer_alt2]
    cohort2 = list(recombination.cross(*parents))

    assert len(cohort2) == 15
    _test_cohort(cohort2, parents, False)
