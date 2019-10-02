import stk
import os
from os.path import join


odir = 'plots'
if not os.path.exists(odir):
    os.mkdir(odir)


def test_progress_plotter(progress):
    plotter1 = stk.ProgressPlotter(
        filename=join(odir, 'fitness'),
        property_fn=lambda mol: mol.fitness,
        y_label='Fitness',
    )
    plotter1.plot(progress)

    plotter2 = stk.ProgressPlotter(
        filename=join(odir, 'atoms'),
        property_fn=lambda mol: len(mol.atoms),
        y_label='Atom Number',
    )
    plotter2.plot(progress)

    # Make a plot where all fitness values are filtered out.
    plotter3 = stk.ProgressPlotter(
        filename=join(odir, 'no_fitness'),
        property_fn=lambda mol: mol.fitness,
        y_label='No Fitness',
        filter=lambda mol: False,
    )
    plotter3.plot(progress)

    # Make a plot where some fitness values are filtered out.
    # Ie by plotting only odd fitness values.
    plotter4 = stk.ProgressPlotter(
        filename=join(odir, 'odd_fitness'),
        property_fn=lambda mol: mol.fitness,
        y_label='Odd Fitness',
        filter=lambda mol: mol.fitness % 1 == 0,
    )
    plotter4.plot(progress)


def test_selection_plotter(progress):
    roulette = stk.Roulette(num_batches=10)
    stk.SelectionPlotter(join('plots', 'roulette_counter'), roulette)
    list(roulette.select(progress.subpopulations[0]))
    list(roulette.select(progress.subpopulations[0]))
