import pytest
from collections import Counter
import numpy as np
from types import SimpleNamespace
import os
from os.path import join
import stk

odir = 'population_tests_output'
if not os.path.exists(odir):
    os.mkdir(odir)


class Mol:
    def __init__(self, x):
        self.x = x

    def same(self, other):
        return other.x == self.x


def test_init():
    """
    Tests the __init__ method of the stk.Population class.

    """

    def NewMol(): return stk.MacroMolecule.__new__(stk.MacroMolecule)

    pop = stk.Population(

              NewMol(), NewMol(),
              stk.Population(NewMol(), NewMol()),
              stk.Population(),
              stk.Population(stk.Population(NewMol(), NewMol())),
              NewMol()
    )

    assert len(pop) == 7
    assert len(pop.members) == 3
    assert len(pop.populations) == 3


def test_add_members_duplicates(generate_population):
    """
    Members in population added to `members` of the other.

    Duplicate additions should be allowed.

    """

    # Create a population to be added and one to be added to.
    receiver = generate_population(cache=True)
    supplier = generate_population(cache=True)

    # Add all mols in `supplier` to `receiver`.
    receiver.add_members(supplier, duplicates=True)

    # Mols in `supplier` should now be held in `members` attribute
    # of `receiver`. A new ``stk.Population`` instance must be created
    # for these tests because `members` is a ``list`` which have a
    # different definition of `__contains`__ to ``stk.Population``.
    # Identity needs to be compared but with lists equality
    # (via ``__eq__``) is compared.
    receiver_members = stk.Population(*receiver.members)
    assert all(mol in receiver_members for mol in supplier)

    # The reverse should not be true.
    supplier_members = stk.Population(*supplier.members)
    assert not all(mol in supplier_members for mol in receiver)

    # Count the frequency of each cage` in `receiver.members`.
    count = Counter(receiver.members)

    # Some should be present twice.
    assert 2 in count.values()

    # No other frequency should be present.
    assert all(freq == 1 or freq == 2 for freq in count.values())


def test_add_members_no_duplicates(generate_population):
    """
    Members in population added to `members` of another.

    Duplicate additions not allowed.

    """

    # Generate an initial populaiton.
    receiver = generate_population(cache=True)

    # Note size of receiver population.
    receiver_size = len(receiver)
    supplier_same = generate_population(cache=True)
    # Add supplier to the receiver. None of the suppliers cages should
    # be added and therefore the len of supplier should be the same as
    # at the start.
    receiver.add_members(supplier_same)
    assert receiver_size == len(receiver)

    supplier_different = generate_population(cache=True, offset=True)
    # Add `supplier_different` to `receiver`. All of the cages should
    # be addable as none should be duplicates. The size of the
    # `receiver` population should increase by the size of the
    # `supplier_different` population.
    receiver.add_members(supplier_different)
    assert receiver_size + len(supplier_different) == len(receiver)


def test_add_subpopulation(generate_population):
    """
    Add a population as a subpopulation to another.

    """

    pop1 = generate_population(cache=True)
    pop2 = generate_population(cache=True)
    pop1.add_subpopulation(pop2)
    assert pop2 not in pop1.populations
    assert all(x in pop1 for x in pop2)


def test_has_structure():

    a1, a2 = Mol(1), Mol(1)
    b1, b2 = Mol(2), Mol(2)

    pop3 = stk.Population()
    pop3.members.append(a1)
    pop3.populations.append(stk.Population())
    pop3.populations[0].members.append(b1)

    assert a1 in pop3
    assert a2 not in pop3
    assert pop3.has_structure(a2)
    assert b1 in pop3
    assert b2 not in pop3
    assert pop3.has_structure(b2)


def test_load(generate_population):
    # Other tests except the cache to be turned off. Make sure it is
    # off after this test.
    try:
        pop = generate_population(cache=True)
        pname = join(odir, 'pop.json')
        pop.dump(pname)

        stk.CACHE_SETTINGS['ON'] = True
        p0 = stk.Population.load(pname, stk.Molecule.from_dict)

        stk.CACHE_SETTINGS['ON'] = False
        p1 = stk.Population.load(pname, stk.Molecule.from_dict)

        for m in p1:
            assert m not in p0

        stk.CACHE_SETTINGS['ON'] = True
        p2 = stk.Population.load(pname, stk.Molecule.from_dict)

        for m in p2:
            assert m in p0

    except Exception:
        raise

    finally:
        stk.CACHE_SETTINGS['ON'] = False


def test_all_members(generate_population, mol):
    """
    Check that all members, direct and in subpopulations, are returned.

    """

    pop = generate_population(cache=True)
    all_members = stk.Population(*(mol for mol in pop.all_members()))

    member_count = 0
    stack = [pop]
    while stack:
        current_pop = stack.pop()
        member_count += len(current_pop.members)
        stack.extend(current_pop.populations)

    assert len(all_members) == member_count


def test_max(generate_population):
    pop = generate_population(cache=False)

    for i, mem in enumerate(pop):
        mem.fitness = i
        mem.unscaled_fitness = [i, 2*i, 3*i, 4*i]

    maxf = pop.max(lambda x: x.fitness)
    maxuf = pop.max(lambda x: x.unscaled_fitness)
    m = np.matrix([x.unscaled_fitness for x in pop])

    assert np.max([x.fitness for x in pop]) == maxf
    assert np.allclose(np.max(m, axis=0), maxuf, atol=1e-8)


def test_mean(generate_population):
    pop = generate_population(cache=False)

    for i, mem in enumerate(pop):
        mem.fitness = i
        mem.unscaled_fitness = [i, 2*i, 3*i, 4*i]

    avgf = pop.mean(lambda x: x.fitness)
    avguf = pop.mean(lambda x: x.unscaled_fitness)
    m = np.matrix([x.unscaled_fitness for x in pop])

    assert np.mean([x.fitness for x in pop]) == avgf
    assert np.allclose(np.mean(m, axis=0), avguf, atol=1e-8)


def test_min(generate_population):
    pop = generate_population(cache=False)

    for i, mem in enumerate(pop):
        mem.fitness = i
        mem.unscaled_fitness = [i, 2*i, 3*i, 4*i]

    minf = pop.min(lambda x: x.fitness)
    minuf = pop.min(lambda x: x.unscaled_fitness)
    m = np.matrix([x.unscaled_fitness for x in pop])

    assert np.min([x.fitness for x in pop]) == minf
    assert np.allclose(np.min(m, axis=0), minuf, atol=1e-8)


def test_remove_duplicates_between_subpops(generate_population):
    """
    Ensure that duplicates are correctly removed from a population.

    """

    subpop1 = generate_population(cache=True)
    subpop2 = generate_population(cache=True)
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


def test_remove_duplicates_not_between_subpops(generate_population):
    """
    Ensures duplicates are removed from within subpopulations only.

    """

    # Create a population from two identical subpopulations.
    subpop1 = generate_population(cache=True)
    subpop2 = generate_population(cache=True)
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


def remove_members(generate_population):
    pop = generate_population(cache=True)
    og_length = len(pop)

    for x in range(5):
        pop[x].remove_me = True

    pop.remove_members(lambda x: hasattr('remove_me'))
    assert len(pop) == og_length - 5


def test_getitem(generate_population):
    """
    Test that the '[]' operator is working.

    """

    # [:] should flatten the population, all members from
    # subpopulations should be transferred to the members attribute.
    pop = generate_population(cache=True)
    flat_pop = pop[:]

    # Change the fitness of one of the members. Fitness should not
    # used in the test for ``in`` and therefore the fact that it is
    # differnt should not matter.
    flat_pop[1].fitness = 555

    assert all(mol in flat_pop.members for mol in pop)
    # Verify lack of subpopulations.
    assert not flat_pop.populations
    # An integer index should return a ``MacroMolecule`` instance.
    assert isinstance(pop[5], stk.MacroMolecule)
    # Non integer slice indices are not supported
    with pytest.raises(TypeError):
        pop[5.5]


def test_sub(generate_population):
    """
    Exclude members of one population from another.

    """

    subtractee = generate_population(cache=True)
    subtractor_same = generate_population(cache=True)
    subtractor_different = generate_population(cache=True, offset=True)

    # Removing mols not present in a population should return the
    # same population.
    result_pop1 = subtractee - subtractor_different
    assert all(mol in subtractee for mol in result_pop1)
    assert len(result_pop1) == len(subtractee)

    # Removing mols should also return a flat population, even if
    # none were actually removed.
    assert not result_pop1.populations

    # Removing mols present in a population should get rid of them.
    result_pop2 = subtractee - subtractor_same
    assert len(result_pop2) == 0


def test_add(generate_population):
    """
    Create a new population from two others.

    The added populations should have their internal structure
    presevered. This means that the way their subpopulations are
    structured is not changed.

    """

    addee = generate_population(cache=True)
    adder = generate_population(cache=True)
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


def test_contains(generate_population):
    """
    Ensure the `in` operator works.

    """

    # Make a population from some mols and initialize.
    mols = [mol for mol in generate_population(cache=True)]

    pop = stk.Population(*mols[:-1])

    # Check that a cage that should not be in it is not.
    assert mols[-1] not in pop
    # Check that a cage that should be in it is.
    assert mols[3] in pop

    # Ensure that the cage is found even if it is in a subpopulation.
    subpop_cages = generate_population(cache=True, offset=True)
    pop.add_subpopulation(subpop_cages)
    assert subpop_cages[2] in pop
