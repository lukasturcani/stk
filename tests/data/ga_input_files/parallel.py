# #####################################################################
# Imports.
# #####################################################################

import stk
import rdkit.Chem.AllChem as rdkit

# #####################################################################
# Initial population.
# #####################################################################

carbon = 'C'
building_blocks = [
    stk.StructUnit2.smiles_init(f'[Br]{carbon*i}[Br]', ['bromine'])
    for i in range(100)
]

topologies = [
    stk.Linear('A', [0], 3),
    stk.Linear('A', [0], 6),
    stk.Linear('A', [0], 12)
]

population = stk.GAPopulation.init_random(stk.Polymer,
                                          [building_blocks],
                                          topologies,
                                          25)

# #####################################################################
# Selector for selecting the next generation.
# #####################################################################

generation_selector = stk.SelectorSequence(
    stk.Fittest(num=3, duplicates=False),
    stk.Roulette(num=22, duplicates=False)
)

# #####################################################################
# Selector for selecting parents.
# #####################################################################

crossover_selector = stk.AboveAverage(num=5, batch_size=2)

# #####################################################################
# Selector for selecting molecules for mutation.
# #####################################################################

mutation_selector = stk.SelectorFunnel(
    stk.AboveAverage(num=10, duplicates=False),
    stk.Roulette(num=5)
)

# #####################################################################
# Crosser.
# #####################################################################

crosser = stk.Jumble(num_offspring_building_blocks=3)

# #####################################################################
# Mutator.
# #####################################################################

mutator = stk.RandomMutation(
    stk.RandomTopology(topologies),
    stk.RandomBuildingBlock(building_blocks, lambda mol: True),
    stk.SimilarBuildingBlock(building_blocks, lambda mol: True, False)
)

# #####################################################################
# Optimizer.
# #####################################################################

optimizer = stk.MMFF()

# #####################################################################
# Fitness calculator.
# #####################################################################


def num_atoms(mol, conformer):
    n_atoms = mol.mol.GetNumAtoms()
    # Save the number of atoms in an attribute for later plotting.
    mol.num_atoms = n_atoms
    return n_atoms


fitness_calculator = stk.PropertyVector(num_atoms)

# #####################################################################
# Fitness normalizer.
# #####################################################################

# The PropertyVector fitness calculator will set the fitness as
# [n_atoms] use the Sum() fitness normalizer to convert the fitness to
# just n_atoms^2. The squre is because we use the Power normalizer.
fitness_normalizer = stk.NormalizerSequence(
    stk.Power(2),
    stk.Sum()
)

# #####################################################################
# Exit condition.
# #####################################################################

exiter = stk.NumGenerations(25)

# #####################################################################
# Make plotters.
# #####################################################################

plotters = [
    stk.ProgressPlotter(filename='fitness_plot',
                        attr='fitness',
                        y_label='Fitness',
                        default=1e-4),
    stk.ProgressPlotter(filename='atom_number_plot',
                        attr='num_atoms',
                        y_label='Number of Atoms',
                        default=0)
]

stk.SelectionPlotter(filename='generational_selection',
                     selector=generation_selector)
stk.SelectionPlotter(filename='crossover_selection',
                     selector=crossover_selector)
stk.SelectionPlotter(filename='mutation_selection',
                     selector=mutation_selector)
