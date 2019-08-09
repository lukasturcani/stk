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
        y_label='Fitness'
    )
    plotter1.plot(progress)

    plotter2 = stk.ProgressPlotter(
        filename=join(odir, 'atoms'),
        property_fn=lambda mol: len(mol.atoms),
        y_label='Atom Number'
    )
    plotter2.plot(progress)


def test_selection_plotter(progress):
    roulette = stk.Roulette(num_batches=10)
    stk.SelectionPlotter(join('plots', 'roulette_counter'), roulette)
    list(roulette.select(progress.subpopulations[0]))
    list(roulette.select(progress.subpopulations[0]))
