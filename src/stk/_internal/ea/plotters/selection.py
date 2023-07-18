import pathlib
import typing
from collections import Counter
from collections.abc import Callable, Iterator, Set
from functools import wraps
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from stk._internal.ea.molecule_record import MoleculeRecord
from stk._internal.ea.selection.batch import Batch, BatchKey
from stk._internal.ea.selection.selectors.selector import Selector
from stk._internal.key_makers.inchi_key import InchiKey

plt.switch_backend("agg")

T = typing.TypeVar("T", bound=MoleculeRecord)


class SelectMethod(typing.Protocol[T]):
    def __call__(
        self,
        population: dict[T, float],
        included_batches: Set[BatchKey] | None = None,
        excluded_batches: Set[BatchKey] | None = None,
    ) -> Iterator[Batch[T]]:
        pass


class SelectionPlotter(typing.Generic[T]):
    """
    Plots which molecule records a :class:`.Selector` selects.

    Examples:
        *Plotting Which Molecule Records Got Selected*

        .. testcode:: plotting-which-molecule-records-got-selected

            import stk

            # Make a selector.
            roulette = stk.Roulette(num_batches=10)

            # Make a population.
            population = {
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
                ): i
                for i in range(100)
            }

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
        filename: pathlib.Path | str,
        selector: Selector[T],
        x_label: str = "Molecule: InChIKey - Fitness Value",
        record_label: Callable[
            [MoleculeRecord[Any], float], str
        ] = lambda record, fitness_value: (
            f"{InchiKey().get_key(record.get_molecule())} - {fitness_value}"
        ),
        heat_map_value: Callable[
            [MoleculeRecord[Any], float], float
        ] = lambda record, fitness_value: fitness_value,
        heat_map_label: str = "Fitness",
        order_by: Callable[
            [MoleculeRecord[Any], float], float
        ] = lambda record, fitness_value: fitness_value,
    ) -> None:
        """
        Parameters:

            filename:
                The basename of the files. This means it should not include
                file extensions.

            selector:
                The :class:`.Selector` whose selection of molecule records
                is plotted.

            x_label:
                The label use for the x axis.

            record_label:
                A function which takes a :class:`.MoleculeRecord` and its
                fitness value and returns a string, which is the label used
                for the :class:`.MoleculeRecord` on the plot.

            heat_map_value:
                A function which takes a :class:`.MoleculeRecord` and its
                fitness value and returns a value. The value is used for
                coloring the heat map used in the plot.

            heat_map_label:
                The label used for the heat map key.

            order_by:
                A function which takes a :class:`.MoleculeRecord` and its
                fitness value and returns a value. The value is used to
                sort the plotted records along the x-axis in descending
                order.
        """
        self._plots = 0
        self._filename = filename
        self._x_label = x_label
        self._record_label = record_label
        self._order_by = order_by
        self._heat_map_value = heat_map_value
        self._heat_map_label = heat_map_label
        selector.select = self._update_counter(selector.select)  # type: ignore

    def _update_counter(self, select: SelectMethod[T]) -> SelectMethod[T]:
        """
        Decorate :meth:`.Selector.select`.

        This is a decorator which makes sure that every time
        :meth:`.Selector.select` selects a :class:`.MoleculeRecord` a
        counter keeping track of selected records is updated.

        Parameters:
            select:
                The :meth:`Selector.select` method to decorate.

        Returns:
            The decorated :meth:`.Selector.select` method.

        """

        @wraps(select)
        def inner(
            population: dict[T, float],
            included_batches: Set[BatchKey] | None = None,
            excluded_batches: Set[BatchKey] | None = None,
        ) -> Iterator[Batch[T]]:
            counter = Counter({record: 0 for record in population})
            for selected in select(
                population=population,
                included_batches=included_batches,
                excluded_batches=excluded_batches,
            ):
                counter.update(selected)
                yield selected
            self._plot(population, counter)

        return inner

    def _plot(self, population: dict[T, float], counter: Counter) -> None:
        """
        Plot a selection counter.

        Parameters:
            counter:
                A counter specifying which records were selected and how
                many times.
        """
        self._plots += 1
        sns.set(style="darkgrid")
        data = []
        for record, selection_count in counter.items():
            label = self._record_label(record, population[record])
            data.append(
                pd.DataFrame(
                    data={
                        self._x_label: label,
                        "Number of Times Selected": selection_count,
                        "order": self._order_by(record, population[record]),
                        "heat_map": self._heat_map_value(
                            record,
                            population[record],
                        ),
                    },
                    index=[self._x_label],
                )
            )
        df = pd.concat(data, ignore_index=True)

        df = df.sort_values(
            ["Number of Times Selected", "order"],
            ascending=[False, False],
        )
        norm = plt.Normalize(df["heat_map"].min(), df["heat_map"].max())
        sm = plt.cm.ScalarMappable(cmap="magma_r", norm=norm)
        sm.set_array([])

        df.to_csv(f"{self._filename}_{self._plots}.csv")
        fig, ax = plt.subplots(figsize=(11.7, 8.28))
        sns.scatterplot(
            x="Number of Times Selected",
            y=self._x_label,
            hue="heat_map",
            palette="magma_r",
            data=df,
            s=[200 for _ in range(len(counter.keys()))],
            ax=ax,
        )
        ax.get_legend().remove()
        # https://tinyurl.com/2p9drmkh
        plt.rcParams["axes.grid"] = False
        ax.figure.colorbar(sm, ax=ax).set_label(self._heat_map_label)
        plt.tight_layout()
        fig.savefig(f"{self._filename}_{self._plots}.png", dpi=fig.dpi)
        plt.close("all")
