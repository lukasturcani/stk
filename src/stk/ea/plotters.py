"""
Plotting
========

#. :class:`.ProgressPlotter`
#. :class:`.SelectionPlotter`

Plotting is done by :class:`.Plotter` objects. Plotters are calculators
which produce graphs. They generally do this
by decorating other calculators to hook into them, collect data from
them and then plot the results. To see examples of how to use plotters
work, look at the documentation of the individual plotters. For example
:class:`SelectionPlotter` and :class:`ProgressPlotter`.


.. _`adding plotters`:

Making New Plotters
-------------------

New plotters should inherit :class:`Plotter` and define a
:meth:`~Plotter.plot` method. Note that plotters do no necessarily have
to make a plot when :meth:`~Plotter.plot` is made. For example,
:class:`SelectionPlotter` plots every time the
:class:`.Selector` it is used with becomes exhausted and it's
:meth:`~SelectionPlotter.plot` does nothing.

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
    Plots how a property changes during a GA run.

    The produced plot will show the GA generations on the x axis and
    the min, mean and max values of an attribute on the y axis.

    Examples
    --------
    Plot how fitness value changes with GA generations

    .. code-block:: python

        import stk

        # progress has subpopulations where each subpopulation is a
        # generation of a GA.
        progress = stk.Population(...)

        # Make the plotter which plots the fitness change across
        # generations.
        plotter = stk.ProgressPlotter(
            filename='fitness_plot',
            property_fn=lambda mol: mol.fitness,
            y_label='Fitness'
        )
        plotter.plot(progress)

    Plot how the number of atoms changes with GA generations

    .. code-block:: python

        plotter = ProgressPlotter(
            filename='atom_plot',
            property_fn=lambda mol: len(mol.atoms),
            y_label='Number of Atoms'
        )
        plotter.plot(progress)

    """

    def __init__(
        self,
        filename,
        property_fn,
        y_label,
        progress_fn=None,
        filter=lambda progress, mol: True,
    ):
        """
        Initialize a :class:`ProgressPlotter` instance.

        Parameters
        ----------
        filename : :class:`str`
            The basename of the files. This means it should not include
            file extensions.

        property_fn : :class:`callable`
            A :class:`callable` which takes a :class:`.EAPopulation`
            and a :class:`.Molecule` object and returns a property
            value of that molecule, which is used for the plot.
            The :class:`callable` must return a valid value for each
            :class:`.Molecule` in the population.

        y_label : :class:`str`
            The y label for the produced graph.

        progress_fn : :class:`callable`, optional
            Takes the population passed to :meth:`plot` and excutes a
            computation on it. This may be useful if you want to
            apply a normalization to the fitness values in the
            progress population, for example.

        filter : :class:`callable`, optional
            Takes an :class:`.EAPopulation` and a :class:`.Molecule` as
            input and returns ``True`` or ``False``. Only molecules
            which return ``True`` will be plotted. Default is for all
            molecules to be plotted.

        """

        self._filename = filename
        self._property_fn = property_fn
        self._y_label = y_label
        self._filter = filter
        self._progress_fn = progress_fn

    def plot(self, progress):
        """
        Plot a progress plot.

        Parameters
        ----------
        progress : :class:`.Population`
            A :class:`.Population` where each generation of the GA is
            a subpopulation.

        Returns
        -------
        None : :class:`NoneType`

        """

        if self._progress_fn is not None:
            self._progress_fn(progress)

        def filter_fn(mol):
            return self._filter(progress, mol)

        def property_fn(mol):
            return self._property_fn(progress, mol)

        sns.set(style='darkgrid')
        df = pd.DataFrame()
        for i, subpop in enumerate(progress.subpopulations, 1):
            filtered = filter(filter_fn, subpop)
            subpop_vals = list(map(property_fn, filtered))

            # If there are no values after filtering, don't plot
            # anything for the generation.
            if not subpop_vals:
                continue

            data = [
                {
                    'Generation': i,
                    self._y_label: max(subpop_vals),
                    'Type': 'Max'
                },
                {
                    'Generation': i,
                    self._y_label: min(subpop_vals),
                    'Type': 'Min'
                },
                {
                    'Generation': i,
                    self._y_label: np.mean(subpop_vals),
                    'Type': 'Mean'
                }
            ]

            df = df.append(data, ignore_index=True)

        # Save the plot data.
        df.to_csv(f'{self._filename}.csv')

        fig = plt.figure(figsize=[8, 4.5])

        palette = sns.color_palette('deep')
        # It's possible that all values were filtered out, and trying
        # to plot an empty dataframe would raise an exception.
        if len(df) != 0:
            sns.scatterplot(
                x='Generation',
                y=self._y_label,
                hue='Type',
                palette={
                    'Max': palette[3],
                    'Min': palette[0],
                    'Mean': palette[2]
                },
                data=df
            )
        # Set the length of the axes to account for all generations,
        # as its possible the first or last ones were not included
        # due to being filtered out.
        plt.xlim(0, len(progress.subpopulations)+1)

        plt.legend(bbox_to_anchor=(1.15, 1), prop={'size': 9})
        plt.tight_layout()
        fig.savefig(f'{self._filename}.png', dpi=500)
        plt.close('all')


class SelectionPlotter(Plotter):
    """
    Plots which molecules a :class:`.Selector` selects.

    Examples
    --------
    .. code-block:: python

        import stk

        # Make a population of molecules.
        pop = stk.Population(...)

        # Make a selector.
        roulette = stk.Roulette(num=10)

        # Make a plotter.
        stk.SelectionPlotter('roulette_counter', roulette)

        # Select the molecules.
        selected = list(roulette.select(pop))

        # There should now be a file called "roulette_counter_1.png"
        # which shows a graph of all the selected molecules.

        # Make another population.
        pop2 = stk.Population(...)

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

    def __init__(
        self,
        filename,
        selector,
        x_label='Molecule: name - fitness value',
        molecule_label=lambda population, mol:
            f'{mol} - {population.get_fitness_values()[mol]}',
        heat_map_value=lambda population, mol:
            population.get_fitness_values()[mol],
        heat_map_label='Fitness',
        order_by=lambda population, mol:
            population.get_fitness_values()[mol],
    ):
        """
        Initialize a :class:`SelectionPlotter` instance.

        Parameters
        ----------
        filename : :class:`str`
            The basename of the files. This means it should not include
            file extensions.

        selector : :class:`.Selector`
            The :class:`.Selector` whose selection of molecules is
            plotted.

        x_label : :class:`str`, optional
            The label use for the x axis.

        molecule_label : :class:`callable`, optional
            A :class:`callable` which takes a :class:`.EAPopulation`
            and a :class:`.Molecule`, for each molecule which is to be
            included on the x-axis of the counter plot. It should
            return a string, which is the label used for the
            :class:`.Molecule` on the plot.

        heat_map_value : :class:`callable`, optional
            A :class:`callable`, which takes a :class:`.EAPopulation`
            and a :class:`.Molecule`, for each molecule which is to be
            included on the x-axis, and returns a value. The value is
            used for coloring the heat map used in the plot.

        heat_map_label : :class:`str`, optional
            The label used for the heat map key.

        order_by : :class:`callable`, optional
            A :class:`callable`, which takes a :class:`.EAPopulation`
            and a :class:`.Molecule`, for each molecule which is to be
            included on the x-axis, and returns a value. The value is
            used to sort the plotted molecules along the x-axis in
            descending order.

        """

        self._plots = 0
        self._filename = filename
        self._x_label = x_label
        self._molecule_label = molecule_label
        self._order_by = order_by
        self._heat_map_value = heat_map_value
        self._heat_map_label = heat_map_label
        selector.select = self._update_counter(selector.select)

    def _update_counter(self, select):
        """
        Decorate :meth:`.Selector.select`.

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
        def inner(population, *args, **kwargs):

            counter = Counter({mol: 0 for mol in population})
            for selected in select(population, *args, **kwargs):
                counter.update(selected)
                yield selected
            self._plot(population, counter)

        return inner

    def _plot(self, population, counter):
        """
        Plot a selection counter.

        Parameters
        ----------
        population : :class:`.EAPopulation`
            The population from which molecules were selected.

        counter : :class:`collections.Counter`
            A counter specifying which molecules were selected and how
            many times.

        Returns
        -------
        None : :class:`NoneType`

        """

        def molecule_label(mol):
            return self._molecule_label(population, mol)

        def heat_map_value(mol):
            return self._heat_map_value(population, mol)

        def order_by(mol):
            return self._order_by(population, mol)

        self._plots += 1
        sns.set(style='darkgrid')
        fig = plt.figure()

        df = pd.DataFrame()
        for mol, selection_count in counter.items():
            label = molecule_label(mol)
            data = {
                self._x_label: label,
                'Number of times selected': selection_count,
                'order': order_by(mol),
                'heat_map': heat_map_value(mol)
            }
            df = df.append(data, ignore_index=True)

        df = df.sort_values(
            ['Number of times selected', 'order'],
            ascending=[False, False]
        )
        norm = plt.Normalize(
            df['heat_map'].min(),
            df['heat_map'].max()
        )
        sm = plt.cm.ScalarMappable(cmap='magma_r', norm=norm)
        sm.set_array([])

        df.to_csv(f'{self._filename}_{self._plots}.csv')
        ax = sns.scatterplot(
            x=self._x_label,
            y='Number of times selected',
            hue='heat_map',
            palette='magma_r',
            data=df,
            s=[200 for i in range(len(counter.keys()))]
        )
        ax.get_legend().remove()
        ax.figure.colorbar(sm).set_label(self._heat_map_label)
        plt.xticks(rotation=90)
        plt.tight_layout()
        fig.savefig(f'{self._filename}_{self._plots}.png', dpi=fig.dpi)
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
