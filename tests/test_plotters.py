import stk
import os
from os.path import join


odir = 'plots'
if not os.path.exists(odir):
    os.mkdir(odir)


def test_progress_plotter(progress):
    plotter1 = stk.ProgressPlotter(filename=join(odir, 'fitness'),
                                   attr='fitness',
                                   y_label='Fitness',
                                   default=None)
    plotter1.plot(progress)

    plotter2 = stk.ProgressPlotter(filename=join(odir, 'cavity_size'),
                                   attr='cavity_size',
                                   y_label='Cavity Size',
                                   default=20)
    plotter2.plot(progress)


def test_selection_plotter(progress):
    roulette = stk.Roulette(num=10)
    stk.SelectionPlotter(join('plots', 'roulette_counter'), roulette)
    list(roulette.select(progress.populations[0]))
    list(roulette.select(progress.populations[0]))
