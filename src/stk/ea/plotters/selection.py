"""
Selection Plotter
=================

"""

from collections import Counter
from functools import wraps
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


plt.switch_backend('agg')


class SelectionPlotter:
    """
    Plots which molecules a :class:`.Selector` selects.

    Examples
    --------
    *Plotting Which Molecules Got Selected*

    .. code-block:: python

        import stk

        # Make a selector.
        roulette = stk.Roulette(num=10)

        # Make a plotter. You do not have to assign it to a variable.
        stk.SelectionPlotter('roulette_counter', roulette)

        # Select the molecules.
        selected = tuple(roulette.select(population))

        # There should now be a file called "roulette_counter_1.png"
        # which shows a graph of all the selected molecules.

        # Select molecules again.
        selected2 = tuple(roulette.select(population))

        # There should now be a file called "roulette_counter_2.png"
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
        :meth:`.Selector.select` selects a :class:`.MoleculeRecord` a
        counter keeping track of selected molecules is updated.

        Parameters
        ----------
        select : :class:`callable`
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
