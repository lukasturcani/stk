import pytest
from collections import Counter
import numpy as np
from types import SimpleNamespace
from os.path import join
import copy

from ..molecular import Cage, MacroMolecule, Molecule
from ..ga import Population, GATools

pop = Population.load(join('data', 'population', 'population.json'),
                      Molecule.fromdict)
pop2 = Population.load('data/population/init_cage_isomers.json',
                       Molecule.fromdict)


def generate_population(offset=False):
    """
    Returns a population of subpopulations and direct members.

    """

    values = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

    # Generate a bunch of cages.
    if offset:
        cages = [Cage.testing_init(values[x+1], values[x],
                 SimpleNamespace(a=values[x])) for x in range(0, 22)]
    if not offset:
        cages = [Cage.testing_init(values[x], values[x],
                 SimpleNamespace(a=values[x])) for x in range(0, 22)]

    # Generate a couple of populations to be used as subpopulations.
    sub1 = Population(*cages[0:4])
    sub2 = Population(*cages[4:9])
    sub3 = Population(*cages[9:14])
    sub4 = Population(*cages[14:19])

    # Place subpopulations into one another.
    sub1.populations.append(sub3)
    sub2.populations.append(sub4)

    # Initialize final population of subpopulations and cages.
    return Population(sub1, sub2, *cages[-3:])


pop = generate_population()


def test_init():
    """
    Tests the __init__ method of the Population class.

    """

    # A test where only ``Cage`` and ``Population`` instances are used
    # for initialization as well as a single ``GATools`` instnace, as
    # intended. Should pass.
    Population(Cage.__new__(Cage), Cage.__new__(Cage),
               Population.__new__(Population), Cage.__new__(Cage),
               Population.__new__(Population), GATools.__new__(GATools))

    # Only ``Cage`` and ``Population`` instances are used for
    # initialization, no ``GATools`` instance. Valid, should pass.
    Population(Cage.__new__(Cage), Cage.__new__(Cage),
               Population.__new__(Population), Cage.__new__(Cage),
               Population.__new__(Population))


def test_add_members_duplicates():
    """
    Members in population added to `members` of the other.

    Duplicate additions should be allowed.

    """

    # Create a population to be added and one to be added to.
    receiver = generate_population()
    supplier = generate_population()

    # Add all cages in `supplier` to `receiver`.
    receiver.add_members(supplier, duplicates=True)

    # Cages in `supplier` should now be held in `members` attribute of
    # `receiver`. A new ``Population`` instance must be created for
    # these tests because `members` is a ``list`` which have a
    # different definition of `__contains`__ to ``Population``.
    # Identity needs to be compared by with lists equality
    # (via ``__eq__``) is compared.
    assert all(cage in Population(*receiver.members)
               for cage in supplier)

    # The reverse should not be true.
    assert not all(cage in Population(*supplier.members)
                   for cage in receiver)

    # Count the frequency of each cage` in `receiver.members`.
    count = Counter(receiver.members)

    # Some should be present twice.
    assert 2 in count.values()

    # No other frequency should be present.
    assert all(freq == 1 or freq == 2 for freq in count.values())


def test_add_members_no_duplicates():
    """
    Members in population added to `members` of another.

    Duplicate additions not allowed.

    """

    # Generate an initial populaiton, initialize the cage's `bb`, `lk`
    # and `topology` attributes.
    receiver = generate_population()

    # Note size of receiver population.
    receiver_size = len(receiver)

    # Same as above but with another population. Note the objects are
    # different instances, but their `lk`, `bb` and `topology`
    # attributes are the same. As a result they should compare equal to
    # those in `receiver`. Indeed, due to caching the same instances
    # should be found in both populations.
    supplier_same = generate_population()

    # Add supplier to the receiver. None of the suppliers cages should
    # be added and therefore the len of supplier should be the same as
    # at the start.
    receiver.add_members(supplier_same)

    assert receiver_size == len(receiver)

    # Generate another population. This time the `bb`, `lk` and
    # `topology` of the cages will have different combinations to the
    # receiver population.
    supplier_different = generate_population(offset=True)

    # Add `supplier_different` to `receiver`. All of the cages should
    # be addable as none should be duplicates. The size of the
    # `receiver` population should increase by the size of the
    # `supplier_different` population.
    receiver.add_members(supplier_different)
    assert receiver_size + len(supplier_different) == len(receiver)


def test_add_subpopulation():
    """
    Add a population as a subpopulation to another.

    """

    pop1 = generate_population()
    pop2 = generate_population()
    pop1.add_subpopulation(pop2)
    assert pop2 not in pop1.populations
    assert all(x in pop1 for x in pop2)


def test_has_structure():

    a1, b1 = pop2[:2]
    a2 = copy.deepcopy(a1)
    b2 = copy.deepcopy(b1)

    pop3 = Population()
    pop3.members.append(a1)
    pop3.populations.append(Population())
    pop3.populations[0].members.append(b1)

    assert a1 in pop3
    assert a2 not in pop3
    assert pop3.has_structure(a2)
    assert b1 in pop3
    assert b2 not in pop3
    assert pop3.has_structure(b2)


def test_load():
    og_cache = dict(Cage.cache)

    pname = join('data', 'population', 'pop2.json')
    p = Population.load(pname,
                        lambda x: Molecule.fromdict(x, load_names=False))
    for mem in p:
        assert not mem.name

    Cage.cache = og_cache

    p = Population.load(pname, Molecule.fromdict)
    for mem in p:
        assert mem.name


def test_all_members():
    """
    Check that all members, direct and in subpopulations, are returned.

    """

    # Generate a bunch of cages.
    cages = [Cage.testing_init(x, 'a', SimpleNamespace(a=1)) for x in
             range(0, 22)]

    # Generate a couple of ``Populations`` to be used as
    # subpopulations.
    sub1 = Population(*cages[0:4])
    sub2 = Population(*cages[4:9])
    sub3 = Population(*cages[9:14])
    sub4 = Population(*cages[14:19])

    # Place subpopulations in one another.
    sub1.populations.append(sub3)
    sub2.populations.append(sub4)

    # Initialize main population from subpopulations and cages.
    main = Population(sub1, sub2, *cages[-3:])

    # Place the results of ``all_members`` into a list.
    all_members = Population(*[cage for cage in main.all_members()])

    # Check that each generated cage is in `all_members`. Should pass.
    assert all(cage in all_members for cage in cages)

    # Add a cage to `cages`. Now there should be a cage in `cages`, not
    # present in main. Should fail.
    cages.append(Cage.testing_init('alpha', 'beta',
                                   SimpleNamespace(a='gamma')))

    with pytest.raises(AssertionError):
        assert all(cage in all_members for cage in cages)


def test_max():
    pop = generate_population()

    for i, mem in enumerate(pop):
        mem.fitness = i
        mem.unscaled_fitness = [i, 2*i, 3*i, 4*i]

    maxf = pop.max(lambda x: x.fitness)
    maxuf = pop.max(lambda x: x.unscaled_fitness)
    m = np.matrix([x.unscaled_fitness for x in pop])

    assert np.max([x.fitness for x in pop]) == maxf
    assert np.allclose(np.max(m, axis=0), maxuf, atol=1e-8)


def test_mean():
    pop = generate_population()

    for i, mem in enumerate(pop):
        mem.fitness = i
        mem.unscaled_fitness = [i, 2*i, 3*i, 4*i]

    avgf = pop.mean(lambda x: x.fitness)
    avguf = pop.mean(lambda x: x.unscaled_fitness)
    m = np.matrix([x.unscaled_fitness for x in pop])

    assert np.mean([x.fitness for x in pop]) == avgf
    assert np.allclose(np.mean(m, axis=0), avguf, atol=1e-8)


def test_min():
    pop = generate_population()

    for i, mem in enumerate(pop):
        mem.fitness = i
        mem.unscaled_fitness = [i, 2*i, 3*i, 4*i]

    minf = pop.min(lambda x: x.fitness)
    minuf = pop.min(lambda x: x.unscaled_fitness)
    m = np.matrix([x.unscaled_fitness for x in pop])

    assert np.min([x.fitness for x in pop]) == minf
    assert np.allclose(np.min(m, axis=0), minuf, atol=1e-8)


def test_remove_duplicates_between_subpops():
    """
    Ensure that duplicates are correctly removed from a population.

    """

    subpop1 = generate_population()
    subpop2 = generate_population()
    main = subpop1 + subpop2

    # Show that cages are duplicated.
    main_count = Counter(main)
    assert all(val == 2 for val in main_count.values())

    # Show that duplicates are not in the same subpopulations.
    subpop_count = Counter(main.populations[0])
    assert all(val == 1 for val in subpop_count.values())

    # Remove duplicates regardless of where they are.
    main.remove_duplicates(between_subpops=True)
    main_count = Counter(main)
    assert all(val == 1 for val in main_count.values())

    # Check that internal structure is maintained.
    assert not main.members
    assert main.populations
    assert main.populations[0].populations
    assert main.populations[0].populations[0].populations
    subsubsubpop = main.populations[0].populations[0].populations[0]
    assert not subsubsubpop.populations


def test_remove_duplicates_not_between_subpops():
    """
    Ensures duplicates are removed from within subpopulations only.

    """

    # Create a population from two identical subpopulations.
    subpop1 = generate_population()
    subpop2 = generate_population()
    main = subpop1 + subpop2

    # Verify that duplicates are present.
    count = Counter(main)
    assert all(val == 2 for val in count.values())

    # Removing duplicates should not change the size of the population
    # as all the duplicates are in different subpopulations.
    main_size = len(main)
    main.remove_duplicates(False)
    assert len(main) == main_size

    # Add one of the subpopulations in the `members` attribute of
    # another, while allowing duplicates.
    subpop1.add_members(subpop2, duplicates=True)
    subpop1_size = len(subpop1)

    # Verify duplicates are present in `members` and in general.
    count2 = Counter(subpop1)
    count2_members = Counter(subpop1.members)
    assert all(val == 2 for val in count2.values())
    assert 2 in count2_members.values()

    # Find the number of duplicates in `members`.
    num_dupes = 0
    for x in count2_members.values():
        if x == 2:
            num_dupes += 1

    # Remove only duplicates within the same subpopulations.
    subpop1.remove_duplicates(between_subpops=False)
    # Size should decrease by the number of duplicates in `members`.
    assert len(subpop1) == subpop1_size - num_dupes

    # Check that internal structure is maintained.
    assert subpop1.members
    assert subpop1.populations
    assert subpop1.populations[0].members
    assert subpop1.populations[0].populations
    assert subpop1.populations[0].populations[0].members
    assert not subpop1.populations[0].populations[0].populations


def remove_members():
    pop = generate_population()
    og_length = len(pop)

    for x in range(5):
        pop[x].remove_me = True

    pop.remove_members(lambda x: hasattr('remove_me'))
    assert len(pop) == og_length - 5


def test_getitem():
    """
    Test that the '[]' operator is working.

    """

    # [:] should flatten the population, all members from
    # subpopulations should be transferred to the members attribute.
    pop = generate_population()
    flat_pop = pop[:]

    # Change the fitness of one of the memebers. Fitness should not
    # used in the test for ``in`` and therefore the fact that it is
    # differnt should not matter.
    flat_pop[1].fitness = 555

    assert all(cage in flat_pop.members for cage in pop)
    # Verify lack of subpopulations.
    assert not flat_pop.populations
    # An integer index should return a ``MacroMolecule`` instance.
    assert isinstance(pop[5], MacroMolecule)
    # Non integer/slice indices are not supported
    with pytest.raises(TypeError):
        pop[5.5]


def test_sub():
    """
    Exclude members of one population from another.

    """

    subtractee = generate_population()
    subtractor_same = generate_population()
    subtractor_different = generate_population(offset=True)

    # Removing cages not present in a population should return the
    # same population.
    result_pop1 = subtractee - subtractor_different
    assert all(cage in subtractee for cage in result_pop1)
    assert len(result_pop1) == len(subtractee)

    # Removing cages should also return a flat population, even if
    # none were actually removed.
    assert not result_pop1.populations

    # Removing cages present in a population should get rid of them.
    result_pop2 = subtractee - subtractor_same
    assert len(result_pop2) == 0


def test_add():
    """
    Create a new population from two others.

    The added populations should have their internal structure
    presevered. This means that the way their subpopulations are
    structured is not changed.

    """

    addee = generate_population()
    adder = generate_population()
    result = addee + adder
    assert len(result) == len(addee) + len(adder)

    # Check that internal structure is maintained
    assert not result.members
    assert result.populations
    assert result.populations[0].members
    assert result.populations[0].populations
    assert result.populations[0].populations[0].members
    assert result.populations[0].populations[0].populations
    assert result.populations[0].populations[0].populations[0].members
    subsubsub_pop = result.populations[0].populations[0].populations[0]
    assert not subsubsub_pop.populations


def test_contains():
    """
    Ensure the `in` operator works.

    """

    # Make a population from some cages and initialize.
    cages = [cage for cage in generate_population()]

    pop = Population(*cages[:-1])

    # Check that a cage that should not be in it is not.
    assert cages[-1] not in pop
    # Check that a cage that should be in it is.
    assert cages[3] in pop

    # Ensure that the cage is found even if it is in a subpopulation.
    subpop_cages = generate_population()

    pop.add_subpopulation(Population(*subpop_cages))
    assert subpop_cages[2] in pop
