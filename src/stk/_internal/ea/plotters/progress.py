import pathlib
import typing
from collections.abc import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

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
                initial_population=[
                    stk.MoleculeRecord(
                        topology_graph=stk.polymer.Linear(
                            building_blocks=[
                                stk.BuildingBlock(
                                    smiles='BrCCBr',
                                    functional_groups=stk.BromoFactory(),
                                ),
                            ],
                            repeating_unit='A',
                            num_repeating_units=i,
                        ),
                    )
                    for i in range(2, 22)
                ],
                fitness_calculator=stk.FitnessFunction(
                    fitness_function=lambda record:
                        record.get_molecule().get_num_atoms(),
                ),
                mutator=stk.RandomBuildingBlock(
                    building_blocks=[
                        stk.BuildingBlock('BrC[Si]CCBr', stk.BromoFactory()),
                        stk.BuildingBlock('BrCCCCCCCBr', stk.BromoFactory()),
                    ],
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

            fitness_values = []
            for generation in ea.get_generations(10):
                fitness_values.append(
                    [
                        fitness_value.raw
                        for fitness_value
                        in generation.get_fitness_values().values()
                    ]
                )

            # Make the plotter which plots the fitness change across
            # generations.
            progress = stk.ProgressPlotter(
                property=fitness_values,
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
                initial_population=[
                    stk.MoleculeRecord(
                        topology_graph=stk.polymer.Linear(
                            building_blocks=[
                                stk.BuildingBlock(
                                    smiles='BrCCBr',
                                    functional_groups=stk.BromoFactory(),
                                ),
                            ],
                            repeating_unit='A',
                            num_repeating_units=i,
                        ),
                    )
                    for i in range(2, 22)
                ],
                fitness_calculator=stk.FitnessFunction(
                    fitness_function=lambda record:
                        record.get_molecule().get_num_atoms(),
                ),
                mutator=stk.RandomBuildingBlock(
                    building_blocks=[
                        stk.BuildingBlock('BrC[Si]CCBr', stk.BromoFactory()),
                        stk.BuildingBlock('BrCCCCCCCBr', stk.BromoFactory()),
                    ],
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

            num_atoms = []
            for generation in ea.get_generations(10):
                num_atoms.append(
                    [
                        record.get_molecule().get_num_atoms()
                        for record in generation.get_molecule_records()
                    ]
                )

            # Make the plotter which plots the number of atoms across
            # generations.
            progress = stk.ProgressPlotter(
                property=num_atoms,
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
                initial_population=[
                    stk.MoleculeRecord(
                        topology_graph=stk.polymer.Linear(
                            building_blocks=[
                                stk.BuildingBlock(
                                    smiles='BrCCBr',
                                    functional_groups=stk.BromoFactory(),
                                ),
                            ],
                            repeating_unit='A',
                            num_repeating_units=i,
                        ),
                    )
                    for i in range(2, 22)
                ],
                fitness_calculator=stk.FitnessFunction(
                    fitness_function=lambda record:
                        record.get_molecule().get_num_atoms(),
                ),
                mutator=stk.RandomBuildingBlock(
                    building_blocks=[
                        stk.BuildingBlock('BrC[Si]CCBr', stk.BromoFactory()),
                        stk.BuildingBlock('BrCCCCCCCBr', stk.BromoFactory()),
                    ],
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

            fitness_values = []
            for generation in ea.get_generations(10):
                fitness_values.append(
                    [
                        fitness_value.normalized
                        for fitness_value
                        in generation.get_fitness_values().values()
                        if fitness_value.raw is not None
                    ]
                )

            # Make the plotter which plots the fitness change across
            # generations.
            progress = stk.ProgressPlotter(
                property=fitness_values,
                y_label='Fitness',
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
        property: Iterable[Sequence[float]],
        y_label: str,
    ) -> None:
        """
        Parameters:
            property (list[list[float]]):
                The generations of the EA, which are plotted.
            y_label:
                The y label for the produced graph.
        """
        self._property = tuple(property)
        self._y_label = y_label
        self._filter = filter
        self._plot_data = self._get_plot_data()

    def _get_plot_data(self) -> pd.DataFrame:
        self._num_generations = 0
        data = []
        for id_, property in enumerate(self._property):
            self._num_generations += 1

            # If there are no values after filtering, don't plot
            # anything for the generation.
            if not property:
                continue

            data.append(
                pd.DataFrame(
                    data={
                        "Generation": [id_, id_, id_],
                        self._y_label: [
                            max(property),
                            np.mean(property),
                            min(property),
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

    def write(self, path: pathlib.Path | str, dpi: int = 500) -> typing.Self:
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
