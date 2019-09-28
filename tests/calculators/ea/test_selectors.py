import stk
import numpy as np
import itertools as it


def _test_no_duplicate_batches(selection):
    selected = set()
    for batch in selection:
        assert batch.get_identity_key() not in selected
        selected.add(batch.get_identity_key())


def _test_no_duplicate_mols(selection):
    selected = set()
    for batch in selection:
        assert all(mol not in selected for mol in batch)
        selected.update(batch)


def _test_selection_properties(selection, num_batches, batch_size):
    for i, batch in enumerate(selection, 1):
        assert batch.get_size() == batch_size
    assert i == num_batches


def _test_selection_determinism(selection1, selection2):
    for batch1, batch2 in it.zip_longest(selection1, selection2):
        assert list(batch1) == list(batch2)


def _test_selection(selection, expected, unselected):
    for mol, in selection:
        assert mol in expected
        assert mol not in unselected


def _base_tests(selector_class, generation, use_random_seed):
    num_batches = [3, 5, 10, 15]
    batch_size = [1, 2, 4, 10]
    duplicate_batches = [True, False]
    duplicate_mols = [True, False]

    options = it.product(
        num_batches,
        batch_size,
        duplicate_batches,
        duplicate_mols,
    )

    for batches, size, dup_batches, dup_mols in options:
        if use_random_seed:
            selector = selector_class(
                num_batches=batches,
                batch_size=size,
                duplicate_batches=dup_batches,
                duplicate_mols=dup_mols,
                random_seed=4,
            )
        else:
            selector = selector_class(
                num_batches=batches,
                batch_size=size,
                duplicate_batches=dup_batches,
                duplicate_mols=dup_mols,
            )

        _test_selection_properties(
            selection=selector.select(generation),
            num_batches=batches,
            batch_size=size,
        )

        _test_selection_determinism(
            selection1=selector.select(generation),
            selection2=selector.select(generation),
        )

        if not dup_batches:
            _test_no_duplicate_batches(selector.select(generation))
        if not dup_mols:
            _test_no_duplicate_mols(selector.select(generation))


def test_best(generation):
    _base_tests(stk.Best, generation, False)

    sorted_gen = sorted(
        generation,
        key=lambda m: m.fitness,
        reverse=True
    )
    best = set(sorted_gen[:10])
    rest = set(sorted_gen[10:])
    _test_selection(stk.Best(10), best, rest)


def test_worst(generation):
    _base_tests(stk.Worst, generation, False)

    sorted_gen = sorted(generation, key=lambda m: m.fitness)
    worst = set(sorted_gen[:10])
    rest = set(sorted_gen[10:])
    _test_selection(stk.Worst(10), worst, rest)


def test_roulette(generation):
    _base_tests(stk.Roulette, generation, True)


def test_above_average(generation):
    _base_tests(stk.AboveAverage, generation, False)


def test_selector_sequence(generation):
    elitism = stk.Fittest(num_batches=5)
    above_avg = stk.AboveAverage()
    selected = (
        *(mol for batch in elitism.select(generation) for mol in batch),
        *(mol for batch in above_avg.select(generation) for mol in batch)
    )
    sequence = stk.SelectorSequence(elitism, above_avg)

    for i, (mol, ) in enumerate(sequence.select(generation)):
        assert mol is selected[i]


def test_stochastic_universal_sampling(generation):
    _base_tests(stk.StochasticUniversalSampling, generation, True)


def test_tournament(generation):
    tournament = stk.Tournament(num_batches=2, batch_size=2)
    pop_batches = tournament._batch(generation)
    worst_batch = min(
        pop_batches,
        key=lambda batch: batch[-1]
    )
    batches = enumerate(tournament.select(generation), 1)
    for i, batch in batches:
        assert len(batch) == 2
        # Ensures the worst performing batch is never selected.
        # This is impossible in tournament sampling.
        assert worst_batch is not batch
    assert i == 2

    tournament_no_dupes = stk.Tournament(
        num_batches=5,
        batch_size=1,
        duplicate_batches=False
    )
    tournament_no_dupes_selected = set(
        mol for mol in tournament_no_dupes.select(generation)
    )
    # Assert that no duplicate molecules are in selected.
    assert len(tournament_no_dupes_selected) == 5

    small_pop = generation[:2]
    small_tournament = stk.Tournament(
        num_batches=10,
        batch_size=1,
        duplicate_batches=True
    )
    # Ensure that out of a choice of two, the lowest fitness
    # is never chosen.
    worst_mol = min(
        small_pop,
        key=lambda mol: mol.fitness
    )
    for batch in small_tournament.select(small_pop):
        assert len(batch) == 1
        assert worst_mol not in batch
