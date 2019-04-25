"""
Defines plotters.

Plotters are calculators which produce graphs. They generally do this
by decorating other calculators to hook into them, collect data from
them and then plot the results. To see examples of how to use plotters
look at the documentation of the individual plotters. For example
:class:`SelectionPlotter` and :class:`ProgressPlotter`.

.. _`adding plotters`:

Extending stk: Making new plotters.
-----------------------------------

New plotters should inherit :class:`Plotter` and define a
:meth:`~Plotter.plot` method. Note that plotters do no necessarily have
to make a plot when :meth:`~Plotter.plot` is made. For example,
:class:`SelectionPlotter` plots every time the
:class:`.Selector` it is used with becomes exhausted and it's
:meth:`~SelectionPlotter.plot` does nothing..

"""

from functools import wraps
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from collections import Counter

plt.switch_backend('agg')


class Plotter:
    """
    A base class for plotters.

    """

    def plot(self, progress):
        """
        Plots a graph.

        Parameters
        ----------
        progress : :class:`.Population`
            A :class:`.Population` where each generation of the GA is
            a subpopulation.

        Returns
        -------
        None : :class:`NoneType`

        """
        raise NotImplementedError()


class ProgressPlotter(Plotter):
    """
    Plots how an attribute changes during a GA.

    The produced plot will show the GA generations on the x axis and
    the min, mean and max values of an attribute on the y axis.

    Attributes
    ----------
    filename : :class:`str`
        The basename of the files. This means it should exclude
        file extensions.

    attr : :class:`str`
        The name of the attribute which is plotted. It must be an
        attribute on the :class:`.Molecule` objects made by the GA.

    y_label : :class:`str`
        The y label for the produced graph.

    default : :class:`float`
        When a molecule in the population does not have attribute
        :attr:`attr`, this is the default value to use for it.

    Examples
    --------
    Plot how fitness value changes with GA generations.

    .. code-block:: python

        # progress has subpopulations where each subpopulation is a
        # generation of a GA.
        progress = Population(...)

        # Make the plotter which plots the fitness change across
        # generations.
        plotter = ProgressPlotter(filename='fitness_plot',
                                  attr='fitness',
                                  y_label='Fitness',
                                  default=1e-4)
        plotter.plot(progress)

    Plot how the number of atoms changes with GA generations.

    .. code-block:: python

        # progress has subpopulations where each subpopulation is a
        # generation of a GA.
        progress = Population(...)

        # Make the plotter which plots the number of atoms across
        # generations. Note that the GA must create the num_atoms
        # attribute on the molecules. This can be done in the fitness
        # function, for example.
        plotter = ProgressPlotter(filename='atom_plot',
                                  attr='num_atoms',
                                  y_label='Number of Atoms',
                                  default=0)
        plotter.plot(progress)

    """

    def __init__(self, filename, attr, y_label, default):
        """
        Initializes a :class:`ProgressPlotter` instance.

        Parameters
        ----------
        filename : :class:`str`
            The basename of the files. This means it should not include
            file extensions.

        attr : :class:`str`
            The name of the attribute which is plotted. It must be an
            attribute on the :class:`.Molecule` objects made by the GA.

        y_label : :class:`str`
            The y label for the produced graph.

        default : :class:`float`
            When a molecule in the population does not have attribute
            :attr:`attr`, this is the default value to use for it.

        """

        self.filename = filename
        self.attr = attr
        self.y_label = y_label
        self.default = default

    def plot(self, progress):
        """
        Plots a progress plot.

        Parameters
        ----------
        progress : :class:`.Population`
            A :class:`.Population` where each generation of the GA is
            a subpopulation.

        Returns
        -------
        None : :class:`NoneType`

        """

        sns.set(style='darkgrid')
        df = pd.DataFrame()
        for i, subpop in enumerate(progress.populations, 1):
            # Ignore failed molecules from progress plot.
            pop = [
                getattr(x, self.attr)
                for x in subpop if getattr(x, self.attr) is not None
            ]
            data = [
                {'Generation': i,
                 self.y_label: max(pop) if pop else self.default,
                 'Type': 'Max'},
                {'Generation': i,
                 self.y_label: min(pop) if pop else self.default,
                 'Type': 'Min'},
                {'Generation': i,
                 self.y_label: np.mean(pop) if pop else self.default,
                 'Type': 'Mean'}
            ]

            df = df.append(data, ignore_index=True)

        # Save the plot data.
        df.to_csv(f'{self.filename}.csv')

        fig = plt.figure(figsize=[8, 4.5])

        palette = sns.color_palette('deep')
        sns.scatterplot(
            x='Generation',
            y=self.y_label,
            hue='Type',
            palette={'Max': palette[3],
                     'Min': palette[0],
                     'Mean': palette[2]},
            data=df
        )

        plt.legend(bbox_to_anchor=(1.15, 1), prop={'size': 9})
        plt.tight_layout()
        fig.savefig(f'{self.filename}.png', dpi=500)
        plt.close('all')


class SelectionPlotter(Plotter):
    """
    Plots which molecules a :class:`.Selector` selects.

    Attributes
    ----------
    filename : :class:`str`
        The basename of the files. This means it should not include
        file extensions.

    selector : :class:`.Selector`
        The :class:`.Selector` whose selection of molecules is
        plotted.

    Examples
    --------
    .. code-block:: python

        # Make a population of molecules.
        pop = Population(...)

        # Make a selector.
        roulette = Roulette(num=10)

        # Make a plotter.
        SelectionPlotter('roulette_counter', roulette)

        # Select the molecules.
        selected = list(roulette.select(pop))

        # There should now be a file called "roulette_counter_1.png"
        # which shows a graph of all the selected molecules.

        # Make another population.
        pop2 = Population(...)

        # Select molecules from this other population.
        selected2 = list(roulette.select(pop2))

        # There should now be a file called "roulette_counter_2.png"
        # which shows a graph of all the selected molecules.

        # Select from the original population again.
        selected3 = list(roulette.select(pop))

        # There should now be a file called "roulette_counter_3.png"
        # which shows a graph of all the selected molecules.

        # And so on every time you use "roulette.select()".

    """

    def __init__(self, filename, selector):
        """
        Initializes a :class:`SelectionPlotter` instance.

        Parameters
        ----------
        filename : :class:`str`
            The basename of the files. This means it should not include
            file extensions.

        selector : :class:`.Selector`
            The :class:`.Selector` whose selection of molecules is
            plotted.

        """

        self.plots = 0
        self.filename = filename
        selector.select = self.update_counter(selector.select)

    def update_counter(self, select):
        """
        Decorates :meth:`.Selector.select`.

        This is a decorator which makes sure that every time
        :meth:`.Selector.select` selects a :class:`.Molecule` a
        counter keeping track of selected molecules is updated.

        Parameters
        ----------
        select : :class:`function`
            The :meth:`Selector.select` method to decorate.

        Returns
        -------
        :class:`function`
            The decorated :meth:`.Selector.select` method.

        """

        @wraps(select)
        def inner(population):

            counter = Counter({mol: 0 for mol in population})
            for selected in select(population):
                counter.update(selected)
                yield selected
            self._plot(counter)

        return inner

    def _plot(self, counter):
        """
        Plots a selection counter.

        Parameters
        ----------
        counter : :class:`Counter`
            A counter specifying which molecules were selected and how
            many times.

        Returns
        -------
        None : :class:`NoneType`

        """

        self.plots += 1
        sns.set(style='darkgrid')
        fig = plt.figure()

        df = pd.DataFrame()
        for mol, selection_count in counter.items():
            label = f'{mol.name} - {mol.fitness}'
            data = {
                'Molecule: name - fitness value': label,
                'Number of times selected': selection_count,
                'Fitness': mol.fitness
            }
            df = df.append(data, ignore_index=True)

        df = df.sort_values(['Number of times selected', 'Fitness'],
                            ascending=[False, False])
        norm = plt.Normalize(df['Fitness'].min(), df['Fitness'].max())
        sm = plt.cm.ScalarMappable(cmap='magma_r', norm=norm)
        sm.set_array([])

        ax = sns.scatterplot(
                    x='Molecule: name - fitness value',
                    y='Number of times selected',
                    hue='Fitness',
                    palette='magma_r',
                    data=df,
                    s=[200 for i in range(len(counter.keys()))])
        ax.get_legend().remove()
        ax.figure.colorbar(sm).set_label('Fitness')
        plt.xticks(rotation=90)
        plt.tight_layout()
        fig.savefig(f'{self.filename}_{self.plots}.png', dpi=fig.dpi)
        plt.close('all')

    def plot(self, progress):
        """
        Does nothing.

        This :class:`Plotter` creates a plot each time :attr:`selector`
        is finished selecting molecules.

        Parameters
        ----------
        progress : :class:`.Population`
            A :class:`.Population` where each generation of the GA is
            a subpopulation.

        Returns
        -------
        None : :class:`NoneType`

        """

        return
