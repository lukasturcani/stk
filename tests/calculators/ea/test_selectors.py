import stk
import numpy as np


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


def test_selector_funnel(generation):
    elitism = stk.Fittest(num_batches=5)
    above_avg = stk.AboveAverage()
    funnel = stk.SelectorFunnel(elitism, above_avg)

    elites = set(
        mol for batch in elitism.select(generation) for mol in batch
    )
    above_avg_elites = set(
        mol for batch in above_avg.select([*elites]) for mol in batch
    )
    funnel_selected = set(
        mol for batch in funnel.select(generation) for mol in batch
    )

    for mol in above_avg_elites:
        assert mol in funnel_selected

    for mol in funnel_selected:
        assert mol in above_avg_elites


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
