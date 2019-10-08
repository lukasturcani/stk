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

generation_selector = stk.Sequence(
    stk.Best(
        num_batches=3,
        duplicate_mols=False,
        duplicate_batches=False
    ),
    stk.RemoveBatches(
        remover=stk.Best(
            num_batches=3,
            duplicate_mols=False,
            duplicate_batches=False,
        ),
        selector=stk.Roulette(
            num_batches=22,
            duplicate_mols=False,
            random_seed=random_seed,
        ),
    ),
)

# #####################################################################
# Selector for selecting parents.
# #####################################################################

crossover_selector = stk.AboveAverage(num_batches=5, batch_size=2)

# #####################################################################
# Selector for selecting molecules for mutation.
# #####################################################################

mutation_selector = stk.FilterBatches(
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

mutator = stk.Random(
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
fitness_normalizer = stk.Sequence(
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
        property_fn=lambda progress, mol:
            progress.get_fitness_values()[mol],
        y_label='Fitness',
        progress_fn=lambda progress:
            progress.set_fitness_values_from_calculators(
                fitness_calculator=fitness_calculator,
                fitness_normalizer=fitness_normalizer,
            )
    ),
    stk.ProgressPlotter(
        filename='atom_number_plot',
        property_fn=lambda progress, mol: len(mol.atoms),
        y_label='Number of Atoms',
    )
]

stk.SelectionPlotter(
    filename='generational_selection',
    selector=generation_selector,
    molecule_label=lambda population, mol:
        f'{mol.id} - {population.get_fitness_values()[mol]}',
    x_label='Molecule: id - fitness value'
)
stk.SelectionPlotter(
    filename='crossover_selection',
    selector=crossover_selector,
    molecule_label=lambda population, mol:
        f'{mol.id} - {population.get_fitness_values()[mol]}',
    x_label='Molecule: id - fitness value'
)
stk.SelectionPlotter(
    filename='mutation_selection',
    selector=mutation_selector,
    molecule_label=lambda population, mol:
        f'{mol.id} - {population.get_fitness_values()[mol]}',
    x_label='Molecule: id - fitness value'
)
