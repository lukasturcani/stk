# #####################################################################
# Imports.
# #####################################################################

import stk
import logging

# #####################################################################
# Set logging level.
# #####################################################################

logging_level = logging.DEBUG

# #####################################################################
# Initial population.
# #####################################################################

carbon = 'C'
building_blocks = [
    stk.BuildingBlock(f'[Br]{carbon*i}[Br]', ['bromine'])
    for i in range(2, 27)
]

topology_graphs = [
    stk.polymer.Linear('A', 3),
    stk.polymer.Linear('A', 6),
    stk.polymer.Linear('A', 12)
]

population = stk.EPopulation.init_random(
    building_blocks=[building_blocks],
    topology_graphs=topology_graphs,
    size=25,
    use_cache=True
)

# #####################################################################
# Selector for selecting the next generation.
# #####################################################################

generation_selector = stk.SelectorSequence(
    stk.Fittest(num_batches=3, duplicates=False),
    stk.Roulette(num_batches=22, duplicates=False)
)

# #####################################################################
# Selector for selecting parents.
# #####################################################################

crossover_selector = stk.AboveAverage(num_batches=5, batch_size=2)

# #####################################################################
# Selector for selecting molecules for mutation.
# #####################################################################

mutation_selector = stk.SelectorFunnel(
    stk.AboveAverage(num_batches=10, duplicates=False),
    stk.Roulette(num_batches=5)
)

# #####################################################################
# Crosser.
# #####################################################################

crosser = stk.Jumble(num_offspring_building_blocks=3)

# #####################################################################
# Mutator.
# #####################################################################

mutator = stk.RandomMutation(
    stk.RandomTopologyGraph(topology_graphs),
    stk.RandomBuildingBlock(building_blocks, lambda mol: True),
    stk.SimilarBuildingBlock(building_blocks, lambda mol: True, False)
)

# #####################################################################
# Optimizer.
# #####################################################################

optimizer = stk.NullOptimizer(use_cache=True)

# #####################################################################
# Fitness calculator.
# #####################################################################


def num_atoms(mol):
    return len(mol.atoms)


fitness_calculator = stk.PropertyVector(num_atoms)

# #####################################################################
# Fitness normalizer.
# #####################################################################

# The PropertyVector fitness calculator will set the fitness as
# [n_atoms] use the Sum() fitness normalizer to convert the fitness to
# just n_atoms^0.5. The sqrt is because we use the Power normalizer.
fitness_normalizer = stk.NormalizerSequence(
    stk.Power(0.5),
    stk.Sum()
)

# #####################################################################
# Exit condition.
# #####################################################################

terminator = stk.NumGenerations(25)

# #####################################################################
# Make plotters.
# #####################################################################

plotters = [
    stk.ProgressPlotter(
        filename='fitness_plot',
        property_fn=lambda mol: mol.fitness,
        y_label='Fitness',
    ),
    stk.ProgressPlotter(
        filename='atom_number_plot',
        property_fn=lambda mol: len(mol.atoms),
        y_label='Number of Atoms',
    )
]

stk.SelectionPlotter(
    filename='generational_selection',
    selector=generation_selector
)
stk.SelectionPlotter(
    filename='crossover_selection',
    selector=crossover_selector
)
stk.SelectionPlotter(
    filename='mutation_selection',
    selector=mutation_selector
)
