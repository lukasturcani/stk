import typing
from collections.abc import Callable, Iterable

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from stk._internal.ea.generation import Generation
from stk._internal.ea.molecule_records.molecule import MoleculeRecord

plt.switch_backend("agg")


class ProgressPlotter:
    """
    Plots how a property changes during an EA run.

    The produced plot will show the EA generations on the x axis and
    the min, mean and max values of an attribute on the y axis.

    Examples:

        *Plotting How Fitness Values Change Across Generations*

        .. testcode:: plotting-how-fitness-values-change-across-generations

            import stk

            # Initialize an EA somehow.
            ea = stk.EvolutionaryAlgorithm(
                initial_population=(
                    stk.MoleculeRecord(
                        topology_graph=stk.polymer.Linear(
                            building_blocks=(
                                stk.BuildingBlock(
                                    smiles='BrCCBr',
                                    functional_groups=[stk.BromoFactory()],
                                ),
                            ),
                            repeating_unit='A',
                            num_repeating_units=i,
                        ),
                    )
                    for i in range(2, 22)
                ),
                fitness_calculator=stk.FitnessFunction(
                    fitness_function=lambda record:
                        record.get_molecule().get_num_atoms(),
                ),
                mutator=stk.RandomBuildingBlock(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='BrC[Si]CCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles='BrCCCCCCCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    is_replaceable=lambda building_block: True
                ),
                crosser=stk.GeneticRecombination(
                    get_gene=lambda building_block: 0
                ),
                generation_selector=stk.Best(
                    num_batches=22,
                    duplicate_molecules=False,
                ),
                mutation_selector=stk.Roulette(
                    num_batches=5,
                    random_seed=10,
                ),
                crossover_selector=stk.Roulette(
                    num_batches=5,
                    batch_size=2,
                    random_seed=10,
                ),
                num_processes=1,
            )

            generations = []
            for generation in ea.get_generations(10):
                generations.append(generation)

            # Make the plotter which plots the fitness change across
            # generations.
            progress = stk.ProgressPlotter(
                generations=generations,
                get_property=lambda record: record.get_fitness_value(),
                y_label='Fitness'
            )
            progress.write('fitness_plot.png')

        .. testcode:: plotting-how-fitness-values-change-across-generations
            :hide:

            import os

            assert os.path.exists('fitness_plot.png')
            os.remove('fitness_plot.png')

        *Plotting How a Molecular Property Changes Across Generations*

        As an example, plotting how the number of atoms changes across
        generations

        .. testcode:: plotting-how-a-molecular-property-changes

            import stk

            # Initialize an EA somehow.
            ea = stk.EvolutionaryAlgorithm(
                initial_population=(
                    stk.MoleculeRecord(
                        topology_graph=stk.polymer.Linear(
                            building_blocks=(
                                stk.BuildingBlock(
                                    smiles='BrCCBr',
                                    functional_groups=[stk.BromoFactory()],
                                ),
                            ),
                            repeating_unit='A',
                            num_repeating_units=i,
                        ),
                    )
                    for i in range(2, 22)
                ),
                fitness_calculator=stk.FitnessFunction(
                    fitness_function=lambda record:
                        record.get_molecule().get_num_atoms(),
                ),
                mutator=stk.RandomBuildingBlock(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='BrC[Si]CCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles='BrCCCCCCCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    is_replaceable=lambda building_block: True
                ),
                crosser=stk.GeneticRecombination(
                    get_gene=lambda building_block: 0
                ),
                generation_selector=stk.Best(
                    num_batches=22,
                    duplicate_molecules=False,
                ),
                mutation_selector=stk.Roulette(
                    num_batches=5,
                    random_seed=10,
                ),
                crossover_selector=stk.Roulette(
                    num_batches=5,
                    batch_size=2,
                    random_seed=10,
                ),
                num_processes=1,
            )

            generations = []
            for generation in ea.get_generations(10):
                generations.append(generation)

            # Make the plotter which plots the number of atoms across
            # generations.
            progress = stk.ProgressPlotter(
                generations=generations,
                get_property=lambda record:
                    record.get_molecule().get_num_atoms(),
                y_label='Number of Atoms'
            )
            progress.write('number_of_atoms_plot.png')

        .. testcode:: plotting-how-a-molecular-property-changes
            :hide:

            import os

            assert os.path.exists('number_of_atoms_plot.png')
            os.remove('number_of_atoms_plot.png')

        *Excluding Molecules From the Plot*

        Sometimes, you want to ignore some molecules from the plot you
        make. For example, If the fitness calculation failed on a
        molecule, you not want to include in a plot of fitness.

        .. testcode:: excluding-molecules-from-the-plot

            import stk

            # Initialize an EA somehow.
            ea = stk.EvolutionaryAlgorithm(
                initial_population=(
                    stk.MoleculeRecord(
                        topology_graph=stk.polymer.Linear(
                            building_blocks=(
                                stk.BuildingBlock(
                                    smiles='BrCCBr',
                                    functional_groups=[stk.BromoFactory()],
                                ),
                            ),
                            repeating_unit='A',
                            num_repeating_units=i,
                        ),
                    )
                    for i in range(2, 22)
                ),
                fitness_calculator=stk.FitnessFunction(
                    fitness_function=lambda record:
                        record.get_molecule().get_num_atoms(),
                ),
                mutator=stk.RandomBuildingBlock(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='BrC[Si]CCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                        stk.BuildingBlock(
                            smiles='BrCCCCCCCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    is_replaceable=lambda building_block: True
                ),
                crosser=stk.GeneticRecombination(
                    get_gene=lambda building_block: 0
                ),
                generation_selector=stk.Best(
                    num_batches=22,
                    duplicate_molecules=False,
                ),
                mutation_selector=stk.Roulette(
                    num_batches=5,
                    random_seed=10,
                ),
                crossover_selector=stk.Roulette(
                    num_batches=5,
                    batch_size=2,
                    random_seed=10,
                ),
                num_processes=1,
            )

            generations = []
            for generation in ea.get_generations(10):
                generations.append(generation)

            # Make the plotter which plots the fitness change across
            # generations.
            progress = stk.ProgressPlotter(
                generations=generations,
                get_property=lambda record: record.get_fitness_value(),
                y_label='Fitness',
                # Only plot records whose unnormalized fitness value is not
                # None, which means the fitness calculation did not fail.
                filter=lambda record:
                    record.get_fitness_value(normalized=False) is not None,
            )
            progress.write('fitness_plot.png')

        .. testcode:: excluding-molecules-from-the-plot
            :hide:

            import os

            assert os.path.exists('fitness_plot.png')
            os.remove('fitness_plot.png')

    """

    def __init__(
        self,
        generations: Iterable[Generation],
        get_property: Callable[[MoleculeRecord], float],
        y_label: str,
        filter: Callable[[MoleculeRecord], bool] = lambda record: True,
    ) -> None:
        """
        Parameters:

            generations (list[Generation]):
                The generations of the EA, which are plotted.

            get_property:
                A :class:`callable` which takes a :class:`.MoleculeRecord`
                and returns a property value of that molecule, which is
                used for the plot. The :class:`callable` must return a
                valid value for each
                :class:`.MoleculeRecord` in `generations`.

            y_label:
                The y label for the produced graph.

            filter:
                Takes an :class:`.MoleculeRecord` and returns
                ``True`` or ``False``. Only records which return ``True``
                are included in the plot. By default, all records will be
                plotted.
        """
        self._get_property = get_property
        self._y_label = y_label
        self._filter = filter
        self._plot_data = self._get_plot_data(generations)

    def _get_plot_data(
        self,
        generations: Iterable[Generation],
    ) -> pd.DataFrame:
        self._num_generations = 0
        data = []
        for id_, generation in enumerate(generations):
            self._num_generations += 1

            filtered = filter(
                self._filter,
                generation.get_molecule_records(),
            )
            properties = tuple(map(self._get_property, filtered))

            # If there are no values after filtering, don't plot
            # anything for the generation.
            if not properties:
                continue

            data.append(
                pd.DataFrame(
                    data={
                        "Generation": [id_, id_, id_],
                        self._y_label: [
                            max(properties),
                            np.mean(properties),
                            min(properties),
                        ],
                        "Type": ["Max", "Mean", "Min"],
                    },
                    index=["Generation", "Generation", "Generation"],
                ),
            )
        return pd.concat(data, ignore_index=True)

    def get_plot_data(self) -> pd.DataFrame:
        """
        Get the plot data.

        Returns:
            A data frame holding the plot data.

        """

        return self._plot_data.copy()

    def write(self, path: str, dpi: int = 500) -> typing.Self:
        """
        Write a progress plot to a file.

        Parameters:

            path:
                The path into which the plot is written.

            dpi:
                The dpi of the image.

        Returns:
            ProgressPlotter: The plotter is returned.

        """

        sns.set(style="darkgrid")
        fig = plt.figure(figsize=[8, 4.5])
        palette = sns.color_palette("deep")

        # It's possible that all values were filtered out, and trying
        # to plot an empty dataframe would raise an exception.
        if len(self._plot_data) != 0:
            sns.scatterplot(
                x="Generation",
                y=self._y_label,
                hue="Type",
                palette={
                    "Max": palette[3],
                    "Min": palette[0],
                    "Mean": palette[2],
                },
                data=self._plot_data,
            )
        # Set the length of the axes to account for all generations,
        # as its possible the first or last ones were not included
        # due to being filtered out.
        plt.xlim(0, self._num_generations)

        plt.legend(bbox_to_anchor=(1.15, 1), prop={"size": 9})
        plt.tight_layout()
        fig.savefig(path, dpi=dpi)
        plt.close("all")
        return self
