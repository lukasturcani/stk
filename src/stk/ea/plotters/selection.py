"""
Selection Plotter
=================

"""

from collections import Counter
from functools import wraps

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from stk.molecular import InchiKey

plt.switch_backend('agg')


class SelectionPlotter:
    """
    Plots which molecule records a :class:`.Selector` selects.

    Examples
    --------
    *Plotting Which Molecule Records Got Selected*

    .. testcode:: plotting-which-molecule-records-got-selected

        import stk

        # Make a selector.
        roulette = stk.Roulette(num_batches=10)

        # Make a population.
        population = tuple(
            stk.MoleculeRecord(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='BrCCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    repeating_unit='A',
                    num_repeating_units=2,
                ),
            ).with_fitness_value(i)
            for i in range(100)
        )

        # Make a plotter. You do not have to assign it to a variable.
        stk.SelectionPlotter('roulette_counter', roulette)

        # Select the molecule records.
        selected = tuple(roulette.select(population))

        # There should now be a file called "roulette_counter_1.png"
        # which shows a graph of all the selected records.

        # Select records again.
        selected2 = tuple(roulette.select(population))

        # There should now be a file called "roulette_counter_2.png"
        # which shows a graph of all the selected molecules.

        # And so on every time you use "roulette.select()".

    .. testcode:: plotting-which-molecule-records-got-selected
        :hide:

        import os

        assert os.path.exists('roulette_counter_1.png')
        assert os.path.exists('roulette_counter_1.csv')
        assert os.path.exists('roulette_counter_2.png')
        assert os.path.exists('roulette_counter_2.csv')
        os.remove('roulette_counter_1.png')
        os.remove('roulette_counter_1.csv')
        os.remove('roulette_counter_2.png')
        os.remove('roulette_counter_2.csv')

    """

    def __init__(
        self,
        filename,
        selector,
        x_label='Molecule: InChIKey - Fitness Value',
        record_label=lambda record: (
            f'{InchiKey().get_key(record.get_molecule())} - '
            f'{record.get_fitness_value()}'
        ),
        heat_map_value=lambda record: record.get_fitness_value(),
        heat_map_label='Fitness',
        order_by=lambda record: record.get_fitness_value(),
    ):
        """
        Initialize a :class:`.SelectionPlotter` instance.

        Parameters
        ----------
        filename : :class:`str`
            The basename of the files. This means it should not include
            file extensions.

        selector : :class:`.Selector`
            The :class:`.Selector` whose selection of molecule records
            is plotted.

        x_label : :class:`str`, optional
            The label use for the x axis.

        record_label : :class:`callable`, optional
            A :class:`callable` which takes a :class:`.MoleculeRecord`
            for each record, which is to be
            included on the x-axis of the counter plot. It should
            return a string, which is the label used for the
            :class:`.MoleculeRecord` on the plot.

        heat_map_value : :class:`callable`, optional
            A :class:`callable`, which takes a :class:`.MoleculeRecord`
            for each record, which is to be included on the x-axis, and
            returns a value. The value is used for coloring the heat
            map used in the plot.

        heat_map_label : :class:`str`, optional
            The label used for the heat map key.

        order_by : :class:`callable`, optional
            A :class:`callable`, which takes a :class:`.MoleculeRecord`
            for each record, which is to be included on the x-axis,
            and returns a value. The value is used to sort the plotted
            records along the x-axis in descending order.

        """

        self._plots = 0
        self._filename = filename
        self._x_label = x_label
        self._record_label = record_label
        self._order_by = order_by
        self._heat_map_value = heat_map_value
        self._heat_map_label = heat_map_label
        selector.select = self._update_counter(selector.select)

    def _update_counter(self, select):
        """
        Decorate :meth:`.Selector.select`.

        This is a decorator which makes sure that every time
        :meth:`.Selector.select` selects a :class:`.MoleculeRecord` a
        counter keeping track of selected records is updated.

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

            counter = Counter({record: 0 for record in population})
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
        population : :class:`tuple` of :class:`.MoleculeRecord`
            The population from which molecule records were selected.

        counter : :class:`collections.Counter`
            A counter specifying which records were selected and how
            many times.

        Returns
        -------
        None : :class:`NoneType`

        """

        self._plots += 1
        sns.set(style='darkgrid')
        df = pd.DataFrame()
        for record, selection_count in counter.items():
            label = self._record_label(record)
            data = {
                self._x_label: label,
                'Number of Times Selected': selection_count,
                'order': self._order_by(record),
                'heat_map': self._heat_map_value(record)
            }
            df = df.append(data, ignore_index=True)

        df = df.sort_values(
            ['Number of Times Selected', 'order'],
            ascending=[False, False]
        )
        norm = plt.Normalize(
            df['heat_map'].min(),
            df['heat_map'].max()
        )
        sm = plt.cm.ScalarMappable(cmap='magma_r', norm=norm)
        sm.set_array([])

        df.to_csv(f'{self._filename}_{self._plots}.csv')
        fig, ax = plt.subplots(figsize=(11.7, 8.28))
        sns.scatterplot(
            x='Number of Times Selected',
            y=self._x_label,
            hue='heat_map',
            palette='magma_r',
            data=df,
            s=[200 for i in range(len(counter.keys()))],
            ax=ax,
        )
        ax.get_legend().remove()
        ax.figure.colorbar(sm).set_label(self._heat_map_label)
        plt.tight_layout()
        fig.savefig(f'{self._filename}_{self._plots}.png', dpi=fig.dpi)
        plt.close('all')
