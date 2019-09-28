# #####################################################################
# Imports.
# #####################################################################

import stk
import logging

# #####################################################################
# Pick a random seed.
# #####################################################################

random_seed = 12

# #####################################################################
# Run GA serially.
# #####################################################################

num_processes = 1

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

population = stk.EAPopulation.init_random(
    building_blocks=[building_blocks],
    topology_graphs=topology_graphs,
    size=25,
    use_cache=True,
    random_seed=random_seed
)

# #####################################################################
# Selector for selecting the next generation.
# #####################################################################

generation_selector = stk.SelectorSequence(
    stk.Best(
        num_batches=3,
        duplicate_mols=False,
        duplicate_batches=False
    ),
    stk.Roulette(
        num_batches=22,
        duplicates=False,
        random_seed=random_seed
    )
)

# #####################################################################
# Selector for selecting parents.
# #####################################################################

crossover_selector = stk.AboveAverage(num_batches=5, batch_size=2)

# #####################################################################
# Selector for selecting molecules for mutation.
# #####################################################################

mutation_selector = stk.SelectorFunnel(
    stk.AboveAverage(duplicate_mols=False),
    stk.Roulette(num_batches=5, random_seed=random_seed)
)

# #####################################################################
# Crosser.
# #####################################################################

crosser = stk.Jumble(
    num_offspring_building_blocks=3,
    random_seed=random_seed
)

# #####################################################################
# Mutator.
# #####################################################################

mutator = stk.RandomMutation(
    stk.RandomTopologyGraph(topology_graphs, random_seed=random_seed),
    stk.RandomBuildingBlock(
        building_blocks=building_blocks,
        key=lambda mol: True,
        random_seed=random_seed
    ),
    stk.SimilarBuildingBlock(
        building_blocks=building_blocks,
        key=lambda mol: True,
        duplicate_building_blocks=False,
        random_seed=random_seed
    ),
    random_seed=random_seed
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
    selector=generation_selector,
    molecule_label=lambda mol: f'{mol.id} - {mol.fitness}',
    x_label='Molecule: id - fitness value'
)
stk.SelectionPlotter(
    filename='crossover_selection',
    selector=crossover_selector,
    molecule_label=lambda mol: f'{mol.id} - {mol.fitness}',
    x_label='Molecule: id - fitness value'
)
stk.SelectionPlotter(
    filename='mutation_selection',
    selector=mutation_selector,
    molecule_label=lambda mol: f'{mol.id} - {mol.fitness}',
    x_label='Molecule: id - fitness value'
)
