import pytest
from collections import Counter
import os
from os.path import join
import stk
import itertools as it

odir = 'population_tests_output'
if not os.path.exists(odir):
    os.mkdir(odir)


def NewMol():
    return stk.ConstructedMolecule.__new__(stk.ConstructedMolecule)


def test_init(amine2, aldehyde2, amine3):
    """
    Tests the __init__ method of the stk.Population class.

    """

    pop = stk.Population(
        NewMol(),
        amine2,
        NewMol(),
        stk.Population(NewMol(), NewMol(), aldehyde2),
        stk.Population(),
        stk.Population(amine3, stk.Population(NewMol(), NewMol())),
        NewMol()
    )

    assert len(pop) == 10
    assert len(pop.direct_members) == 4
    assert len(pop.subpopulations) == 3


def test_init_all():

    amines = [
        stk.BuildingBlock('NCCCN', ['amine']),
        stk.BuildingBlock('NCCCCCN', ['amine']),
        stk.BuildingBlock('NCCOCCN', ['amine']),
    ]
    aldehydes = [
        stk.BuildingBlock('O=CCC(C=O)CC=O', ['aldehyde']),
        stk.BuildingBlock('O=CCC(C=O)CC=O', ['aldehyde']),
        stk.BuildingBlock('O=CC(C=O)COCC=O', ['aldehyde']),
    ]
    # A total of 9 cages will be created.
    cages = stk.Population.init_all(
        building_blocks=[amines, aldehydes],
        topology_graphs=[stk.cage.FourPlusSix()]
    )

    assert len(cages) == 9
    cages.remove_duplicates(key=lambda mol: mol.get_identity_key())
    assert len(cages) == 6

    for cage in cages:
        bbs = tuple(cage.building_block_vertices.keys())
        assert len(bbs) == 2
        assert (
            repr(cage.topology_graph) == repr(stk.cage.FourPlusSix())
        )
        amines = tuple(
            bb for bb in bbs
            if bb.func_groups[0].fg_type.name == 'amine'
        )
        assert len(amines) == 1
        aldehydes = tuple(
            bb for bb in bbs
            if bb.func_groups[0].fg_type.name == 'aldehyde'
        )
        assert len(aldehydes) == 1


def test_clone(population):
    clone = population.clone()
    for m1, m2 in it.zip_longest(population, clone):
        assert m1 is m2


def test_init_random():
    amines = [
        stk.BuildingBlock('NCCCN', ['amine']),
        stk.BuildingBlock('NCCCCCN', ['amine']),
        stk.BuildingBlock('NCCOCCN', ['amine']),
    ]
    aldehydes = [
        stk.BuildingBlock('O=CCC(C=O)CC=O', ['aldehyde']),
        stk.BuildingBlock('O=CCC(C=O)CC=O', ['aldehyde']),
        stk.BuildingBlock('O=CC(C=O)COCC=O', ['aldehyde']),
    ]

    cages = stk.Population.init_random(
        building_blocks=[amines, aldehydes],
        topology_graphs=[stk.cage.FourPlusSix()],
        size=5
    )

    assert len(cages) == 5
    assert len(cages.direct_members) == 5
    assert len(cages.subpopulations) == 0

    for cage in cages:
        bbs = tuple(cage.building_block_vertices.keys())
        assert len(bbs) == 2
        assert (
            repr(cage.topology_graph) == repr(stk.cage.FourPlusSix())
        )
        amines = tuple(
            bb for bb in bbs
            if bb.func_groups[0].fg_type.name == 'amine'
        )
        assert len(amines) == 1
        aldehydes = tuple(
            bb for bb in bbs
            if bb.func_groups[0].fg_type.name == 'aldehyde'
        )
        assert len(aldehydes) == 1


def test_init_diverse():
    amines = [
        stk.BuildingBlock('NCCCN', ['amine']),
        stk.BuildingBlock('NCCCCCN', ['amine']),
        stk.BuildingBlock('NCCOCCN', ['amine']),
    ]
    aldehydes = [
        stk.BuildingBlock('O=CCC(C=O)CC=O', ['aldehyde']),
        stk.BuildingBlock('O=CCC(C=O)CC=O', ['aldehyde']),
        stk.BuildingBlock('O=CC(C=O)COCC=O', ['aldehyde']),
    ]

    # A total of 4 cages will be created.
    cages = stk.Population.init_diverse(
        building_blocks=[amines, aldehydes],
        topology_graphs=[stk.cage.FourPlusSix()],
        size=2
    )

    cage1, cage2 = cages
    amine = next(
        bb for bb in cage1.building_block_vertices
        if bb.func_groups[0].fg_type.name == 'amine'
    )
    aldehyde = next(
        bb for bb in cage1.building_block_vertices
        if bb.func_groups[0].fg_type.name == 'aldehyde'
    )
    diff_bb1 = min(
        amines,
        key=lambda m: stk.dice_similarity(m, amine)
    )
    diff_bb2 = min(
        aldehydes,
        key=lambda m: stk.dice_similarity(m, aldehyde)
    )
    cage2_amine = next(
        bb for bb in cage2.building_block_vertices
        if bb.func_groups[0].fg_type.name == 'amine'
    )
    cage2_aldehyde = next(
        bb for bb in cage2.building_block_vertices
        if bb.func_groups[0].fg_type.name == 'aldehyde'
    )
    assert cage2_amine is diff_bb1
    assert cage2_aldehyde is diff_bb2


def test_add_members(tmp_population, population):
    assert len(tmp_population) == len(population)

    tmp_population.add_members(
        molecules=population,
        duplicate_key=lambda m: m.get_identity_key()
    )
    assert len(population) == len(tmp_population)

    tmp_population.add_members(
        molecules=population,
        duplicate_key=None
    )
    assert len(population)*2 == len(tmp_population)


def test_add_subpopulation(population, tmp_population):
    for member in population:
        assert member not in tmp_population

    assert len(population) == len(tmp_population)

    num_direct_members = len(tmp_population.direct_members)
    tmp_population.add_subpopulation(population)
    assert len(tmp_population.direct_members) == num_direct_members
    assert len(population)*2 == len(tmp_population)

    for member in population:
        assert member in tmp_population


def test_get_all_members(population):
    """
    Check that all members, direct and in subpopulations, are returned.

    """

    assert len(list(population)) == len(population)


def test_set_mol_ids(tmp_population):
    for member in tmp_population:
        assert not hasattr(member, 'id')

    tmp_population.set_mol_ids(10)
    ids = set(m.id for m in tmp_population)
    tmp_population.remove_duplicates()
    for i in range(10, 10+len(tmp_population)):
        assert i in ids
    assert len(ids) == len(tmp_population)


def test_dump_and_load(tmp_population):
    path = join(odir, 'population.dump')

    # Add optional attrs to one of the mooecules.
    first = tmp_population[0]
    first.test_attr1 = 'something'
    first.test_attr2 = 12
    first.test_attr3 = ['12', 'something', 21]
    first.test_attr4 = 1111333
    include_attrs = [
        'test_attr1', 'test_attr2', 'test_attr3', 'test_attr5'
    ]

    with pytest.raises(Exception):
        tmp_population.dump(
            path=path,
            include_attrs=include_attrs,
            ignore_missing_attrs=False
        )

    tmp_population.dump(
        path=path,
        include_attrs=include_attrs,
        ignore_missing_attrs=True
    )

    p0 = stk.Population.load(path, use_cache=False)
    # Check that the molecule has the extra attributes.
    assert p0[0].test_attr1 == first.test_attr1
    assert p0[0].test_attr2 == first.test_attr2
    assert p0[0].test_attr3 == first.test_attr3
    assert not hasattr(p0[0], 'test_attr4')

    # Check that other molecules do not have extra attributes.
    for mol in p0[1:]:
        assert all(not hasattr(mol, attr) for attr in include_attrs)

    p1 = stk.Population.load(path, use_cache=True)
    # Check that the molecule has the extra attributes.
    assert p1[0].test_attr1 == first.test_attr1
    assert p1[0].test_attr2 == first.test_attr2
    assert p1[0].test_attr3 == first.test_attr3
    assert not hasattr(p1[0], 'test_attr4')
    assert all(m not in p0 for m in p1)

    # Check that other molecules do not have extra attributes.
    for mol in p0[1:]:
        assert all(not hasattr(mol, attr) for attr in include_attrs)

    p2 = stk.Population.load(path, use_cache=True)
    assert all(m in p2 for m in p1)


def test_optimize(tmp_population):
    optimizer = stk.NullOptimizer(use_cache=True)
    assert not optimizer._cache
    tmp_population.optimize(optimizer)
    assert len(optimizer._cache) == len(set(tmp_population))

    raiser = stk.RaisingCalculator(optimizer, 1)
    with pytest.raises(stk.RaisingCalculatorError):
        tmp_population.optimize(raiser)


def test_remove_duplicates_across_subpopulations(tmp_population):
    tmp_population.remove_duplicates(across_subpopulations=True)
    main = tmp_population.clone() + tmp_population.clone()

    # Show that cages are duplicated.
    assert all(val == 2 for val in Counter(main).values())
    # Show that duplicates are not in the same subpopulations.
    assert all(
        val == 1 for val in Counter(main.subpopulations[0]).values()
    )
    # Remove duplicates regardless of where they are.
    main.remove_duplicates(across_subpopulations=True)
    assert all(val == 1 for val in Counter(main).values())


def test_remove_duplicates_not_across_subpopulations():
    x = stk.BuildingBlock('C')
    y = stk.BuildingBlock('N')
    subpop = stk.Population(x, x, y, stk.Population(y, y))
    pop = subpop.clone() + subpop.clone()
    counts = Counter(pop)
    assert counts[x] == 4
    assert counts[y] == 6

    subpop_counts = Counter(pop.subpopulations[0].direct_members)
    assert subpop_counts[x] == 2
    assert subpop_counts[y] == 1

    pop.remove_duplicates(across_subpopulations=False)
    counts = Counter(pop)
    assert counts[x] == 2
    assert counts[y] == 4

    subpop_counts = Counter(pop.subpopulations[0].direct_members)
    assert subpop_counts[x] == 1
    assert subpop_counts[y] == 1


def remove_members(tmp_population):
    assert any(
        isinstance(m, stk.BuildingBlock) for m in tmp_population
    )
    tmp_population.remove_members(
        key=lambda m: isinstance(m, stk.BuildingBlock)
    )
    assert not any(
        isinstance(m, stk.BuildingBlock) for m in tmp_population
    )


def test_getitem(population):
    """
    Test that the '[]' operator is working.

    """

    # [:] should flatten the population, all members from
    # subpopulations should be transferred to the members attribute.
    flat_pop = population[:]
    assert len(flat_pop.direct_members) == len(population)
    assert population[1] is population.direct_members[1]


def test_len(population):
    assert len(population) == 19


def test_sub(population):
    assert len(population - population.clone()) == 0


def test_add(population):
    """
    Create a new population from two others.

    The added populations should have their internal structure
    presevered. This means that the way their subpopulations are
    structured is not changed.

    """

    c1 = population.clone()
    c2 = population.clone()
    result = c1 + c2
    assert len(result) == 2*len(population)
    assert not result.direct_members
    assert len(result.subpopulations) == 2
    assert result.subpopulations[0] is c1
    assert result.subpopulations[1] is c2


def test_contains(tmp_population, population):
    """
    Ensure the `in` operator works.

    """

    assert all(m not in population for m in tmp_population)
    assert all(m in population for m in population.clone())
