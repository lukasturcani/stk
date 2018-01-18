# Import stuff
from mtk import FunctionData, FourPlusSix

# #################################################################################
# Define variables.
# #################################################################################

bb_db_path = r'path1'
lk_db_path = r'path2'
var1 = FunctionData('func1', bb=bb_db_path, lk=lk_db_path)

# #################################################################################
# Nested variables test.
# #################################################################################

one = '1'
two = one
three = two


# #################################################################################
# Databases of stored molecules to use.
# #################################################################################

databases = [one, '2']

# #################################################################################
# Run mmea serially.
# #################################################################################

processes = 1


# #################################################################################
# Population initialization function.
# #################################################################################

init_func = {'NAME': 'init_random_cages',
             'bb_db': bb_db_path,
             'lk_db': lk_db_path,
             'topologies': [FourPlusSix()]}

# #################################################################################
# Selection function for selecting the next generation.
# #################################################################################

generational_select_func = {'NAME': 'stochastic_sampling',
                            'use_rank': True}

# #################################################################################
# Selection function for selecting parents.
# #################################################################################

crossover_select_func = {'NAME': 'crossover_roulette'}

# #################################################################################
# Selection function for selecting molecules for mutation.
# #################################################################################

mutation_select_func = {'NAME': 'stochastic_sampling',
                        'duplicates': True}

# #################################################################################
# Crossover function.
# #################################################################################

crossover_func1 = {'NAME': 'bb_lk_exchange'}
crossover_funcs = [crossover_func1]

# #################################################################################
# Mutation function 1.
# #################################################################################

mutation_func1 = {'NAME': 'cage_random_bb',
                  'database': bb_db_path}

# #################################################################################
# Mutation function 2.
# #################################################################################

mutation_func2 = {'NAME': 'cage_random_lk',
                  'database': lk_db_path}

mutation_funcs = [mutation_func1, mutation_func2]

# #################################################################################
# When carrying mutations, chance that a given mutation function will be used.
# #################################################################################

mutation_weights = None

# #################################################################################
# Optimization function.
# #################################################################################

opt_func = {'NAME': 'do_not_optimize'}

# #################################################################################
# Fitness function.
# #################################################################################

fitness_func = {'NAME': 'random_fitness_vector'}

# #################################################################################
# Normalization functions.
# #################################################################################

# First shift all energy values so that they are always positive.
normalization_func1 = {'NAME': 'shift_elements',
                       'indices': [-1]}

# Second, make sure that the magnitudes of all fitness parameters
# are comparable.
normalization_func2 = {'NAME': 'magnitudes'}


normalization_funcs = [normalization_func1, normalization_func2]

# #################################################################################
# Number of generations to create.
# #################################################################################

num_generations = 5

# #################################################################################
# Number of mutation operations to perform each generation.
# #################################################################################

num_mutations = 2

# #################################################################################
# Number of crossover operations to perform each generation.
# #################################################################################

num_crossovers = 2

# #################################################################################
# Size of the population.
# #################################################################################

pop_size = 5

comparison_pops = ['pop1', 'pop2']
