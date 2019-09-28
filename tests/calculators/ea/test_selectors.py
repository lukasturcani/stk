import stk
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


def _test_num_batches(selection, num_batches):
    assert sum(1 for _ in selection) == num_batches


def _test_batch_size(selection, batch_size):
    ...


def _test_selection_determinism(selection1, selection2):
    ...


def _test_expected_selection(selection, expected):
    ...


def _test_expected_unselected(selection, unselected):
    ...


def _test_selection_order(selection, expected):
    ...


def test_fittest(generation):
    fittest = stk.Fittest(num_batches=10)
    selected = fittest.select(generation)
    sorted_pop = sorted(
        generation,
        reverse=True,
        key=lambda mol: mol.fitness
    )
    for (mol1, ), mol2 in zip(selected, sorted_pop):
        assert mol1 is mol2

    fittest = stk.Fittest(batch_size=5)
    for batch in fittest.select(generation):
        assert len(batch) == 5


def test_roulette(generation):
    roulette = stk.Roulette(num_batches=5, batch_size=5)
    for batch in roulette.select(generation):
        assert len(batch) == 5


def test_above_average(generation):
    mean = np.mean([mol.fitness for mol in generation])
    above_average = stk.AboveAverage()
    selected = set(
        mol
        for batch in above_average.select(generation)
        for mol in batch
    )

    for mol in generation:
        if mol.fitness >= mean:
            assert mol in selected
        else:
            assert mol not in selected

    above_average = stk.AboveAverage(batch_size=5)
    for batch in above_average.select(generation):
        assert len(batch) == 5


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


def test_stochastic_universal(generation):
    stochastic = stk.StochasticUniversalSampling(num_batches=5, batch_size=5)
    for batch in stochastic.select(generation):
        assert len(batch) == 5
    assert len(list(stochastic.select(generation))) == 5
    stochastic_ranked = stk.StochasticUniversalSampling(
        num_batches=5,
        batch_size=1,
        use_rank=True
    )
    assert len(list(stochastic_ranked.select(generation))) == 5


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
