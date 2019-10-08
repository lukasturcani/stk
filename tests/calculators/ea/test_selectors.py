import stk
import itertools as it
from functools import partial
import numpy as np


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


def _test_expected(selection, expected):
    for mol, in selection:
        assert mol in expected


def _test_unselected(selection, unselected):
    for mol, in selection:
        assert mol not in unselected


def _base_tests(
    selector_class,
    generation,
    num_batches,
    batch_size,
    duplicate_batches,
    duplicate_mols,
    use_random_seed
):
    options = it.product(
        num_batches,
        batch_size,
        duplicate_batches,
        duplicate_mols,
    )

    for batches, size, dup_batches, dup_mols in options:
        if use_random_seed:
            selector = partial(
                selector_class,
                num_batches=batches,
                batch_size=size,
                duplicate_batches=dup_batches,
                duplicate_mols=dup_mols,
                random_seed=4,
            )
        else:
            selector = partial(
                selector_class,
                num_batches=batches,
                batch_size=size,
                duplicate_batches=dup_batches,
                duplicate_mols=dup_mols,
            )

        _test_selection_properties(
            selection=selector().select(generation),
            num_batches=batches,
            batch_size=size,
        )

        _test_selection_determinism(
            selection1=selector().select(generation),
            selection2=selector().select(generation),
        )

        if not dup_batches:
            _test_no_duplicate_batches(selector().select(generation))
        if not dup_mols:
            _test_no_duplicate_mols(selector().select(generation))


def test_best(generation):
    _base_tests(
        selector_class=stk.Best,
        generation=generation,
        num_batches=[3, 5, 10],
        batch_size=[1, 2],
        duplicate_batches=[True, False],
        duplicate_mols=[True, False],
        use_random_seed=False,
    )

    sorted_gen = sorted(
        generation,
        key=lambda m: m.fitness,
        reverse=True
    )
    best = set(sorted_gen[:10])
    rest = set(sorted_gen[10:])
    _test_expected(stk.Best(10).select(generation), best)
    _test_unselected(stk.Best(10).select(generation), rest)


def test_worst(generation):
    _base_tests(
        selector_class=stk.Worst,
        generation=generation,
        num_batches=[3, 5, 10],
        batch_size=[1, 2],
        duplicate_batches=[True, False],
        duplicate_mols=[True, False],
        use_random_seed=False,
    )

    sorted_gen = sorted(generation, key=lambda m: m.fitness)
    worst = set(sorted_gen[:10])
    rest = set(sorted_gen[10:])
    _test_expected(stk.Worst(10).select(generation), worst)
    _test_unselected(stk.Worst(10).select(generation), rest)


def test_roulette(generation):
    _base_tests(
        selector_class=stk.Roulette,
        generation=generation,
        num_batches=[3, 5, 10],
        batch_size=[1, 2],
        duplicate_batches=[True, False],
        duplicate_mols=[True, False],
        use_random_seed=True,
    )


def test_above_average(generation):
    _base_tests(
        selector_class=stk.AboveAverage,
        generation=generation,
        num_batches=[3, 5, 10],
        batch_size=[1, 2],
        duplicate_batches=[True, False],
        duplicate_mols=[True, False],
        use_random_seed=False,
    )
    mean = np.mean([m.fitness for m in generation])
    expected = {m for m in generation if m.fitness > mean}
    unselected = {m for m in generation if m.fitness < mean}
    _test_expected(stk.AboveAverage().select(generation), expected)
    _test_unselected(stk.AboveAverage().select(generation), unselected)


def test_stochastic_universal_sampling(generation):
    _base_tests(
        selector_class=stk.StochasticUniversalSampling,
        generation=generation,
        num_batches=[3, 5, 10],
        batch_size=[1, 2],
        duplicate_batches=[True, False],
        duplicate_mols=[True, False],
        use_random_seed=True,
    )


def test_tournament(generation):
    _base_tests(
        selector_class=stk.Tournament,
        generation=generation,
        num_batches=[3, 5, 10],
        batch_size=[1, 2],
        duplicate_batches=[True, False],
        duplicate_mols=[True, False],
        use_random_seed=True,
    )
    tournament = stk.Tournament(duplicate_mols=False)
    fitness_values = generation.get_fitness_values()
    _test_unselected(
        selection=tournament.select(generation),
        unselected={min(generation, key=lambda m: fitness_values[m])},
    )


def test_remove_batches(generation):
    fitness_values = generation.get_fitness_values()
    sorted_gen = sorted(generation, key=lambda m: fitness_values[m])
    worst = sorted_gen[:10]
    next_worst = sorted_gen[10:20]

    removed = stk.RemoveBatches(
        remover=stk.Worst(10),
        selector=stk.Worst(10),
    )

    _test_expected(removed.select(generation), next_worst)
    _test_unselected(removed.select(generation), worst)


def test_remove_molecules(generation):
    fitness_values = generation.get_fitness_values()
    sorted_gen = sorted(generation, key=lambda m: fitness_values[m])
    worst = sorted_gen[:10]
    next_worst = sorted_gen[10:20]

    removed = stk.RemoveBatches(
        remover=stk.Worst(10),
        selector=stk.Worst(10),
    )

    _test_expected(removed.select(generation), next_worst)
    _test_unselected(removed.select(generation), worst)


def test_filter_batches(generation):
    fitness_values = generation.get_fitness_values()
    sorted_gen = sorted(generation, key=lambda m: fitness_values[m])
    worst = sorted_gen[:10]
    rest = sorted_gen[10:]

    filtered = stk.FilterBatches(
        filter=stk.Worst(10),
        selector=stk.Worst(),
    )

    _test_expected(filtered.select(generation), worst)
    _test_unselected(filtered.select(generation), rest)


def test_filter_molecules(generation):
    fitness_values = generation.get_fitness_values()
    sorted_gen = sorted(generation, key=lambda m: fitness_values[m])
    worst = sorted_gen[:10]
    rest = sorted_gen[10:]

    filtered = stk.FilterMolecules(
        filter=stk.Worst(10),
        selector=stk.Worst(),
    )

    _test_expected(filtered.select(generation), worst)
    _test_unselected(filtered.select(generation), rest)
