# Import stuff.
import os
from os.path import join

##################################################################################
# Define variables.
##################################################################################

# Database variables.
cx1 = '/work/lt912/MMEA/dbs/liverpool_refined'
local = '/Users/marcinmiklitz/liverpool_refined'

path = next(x for x in [cx1, local] if os.path.exists(x))

bb_db = join(path, 'aldehydes_3f')
lk_db = join(path, 'amines_2f')

# Schrodinger variables.
cx1 = '/work/lt912/schrodinger2016-4'
local = '/opt/schrodinger/suites2017-1'

mm_path = next(x for x in [cx1, local] if os.path.exists(x))

##################################################################################
# Population initialization function.
##################################################################################

init_func = {'NAME' : 'init_random_cages',
             'bb_db' : bb_db,
             'lk_db' : lk_db,
             'topologies' : [FourPlusSix()]}

##################################################################################
# Selection function for selecting the next generation.
##################################################################################

generational_select_func =  {'NAME' : 'stochastic_sampling',
                             'use_rank' : True}

##################################################################################
# Selection function for selecting parents.
##################################################################################

parent_select_func =  {'NAME' : 'crossover_roulette'}

##################################################################################
# Selection function for selecting molecules for mutation.
##################################################################################

mutant_select_func = {'NAME' : 'stochastic_sampling',
                      'duplicates' : True}

##################################################################################
# Crossover functions.
##################################################################################

crossover_funcs =  [{'NAME' : 'bb_lk_exchange'}]

##################################################################################
# Mutation function 1.
##################################################################################

mutation_func1 = {'NAME' : 'cage_random_bb',
                 'database' : bb_db}

##################################################################################
# Mutation function 2.
##################################################################################

mutation_func2 = {'NAME' : 'cage_random_lk',
                 'database' : lk_db}

##################################################################################
# Mutation function 3.
##################################################################################

mutation_func3 = {'NAME' : 'cage_similar_bb',
                 'database' : bb_db}

##################################################################################
# Mutation function 4.
##################################################################################

mutation_func4 = {'NAME' : 'cage_similar_lk',
                 'database' : lk_db}

##################################################################################
# Joining mutation functions.
##################################################################################

mutation_funcs = [mutation_func1, mutation_func2, mutation_func3, mutation_func4]

##################################################################################
# When carrying mutations, chance that a given mutation function will be used.
##################################################################################

mutation_weights = [1/4,1/4,1/4,1/4]

##################################################################################
# Optimization function.
##################################################################################

opt_func = {'NAME' : 'macromodel_cage_opt',
            'macromodel_path' : mm_path}

##################################################################################
# Fitness function.
##################################################################################

efunc = FunctionData('macromodel', forcefield=16, macromodel_path=mm_path)
fitness_func = {'NAME' : 'cage',
                'pseudoformation_params' : {'func' : efunc}}

##################################################################################
# Normalization functions.
##################################################################################

norm_func0 = {'NAME' : 'cage',
              'cavity' : 8,
              'window' : 8}

# First shift all energy values so that they are always positive.
norm_func1 = {'NAME' : 'shift_elements',
                      'indices' : [-1]}

# Second, make sure that the magnitudes of all fitness parameters
# are comparable.
norm_func2 = {'NAME' : 'magnitudes'}

# Third, combine the fitness parameter into a single fitness value.
norm_func3 =  {'NAME' : 'combine',
                        'coefficients' : [1,1,1,1],
                        'exponents' : [1,1,1,1]}

# Last, invert the fitness value because all fitness parameters
# are inversely proportional to fitness.
norm_func4 =  {'NAME' : 'invert'}

# Now define the normalization_funcs variable holding all of the
# normalization functions to be used.

normalization_funcs = [norm_func0, norm_func1,
                       norm_func2, norm_func3, norm_func4]

##################################################################################
# Number of generations to create.
##################################################################################

num_generations=2

##################################################################################
# Number of mutation operations to perform each generation.
##################################################################################

num_mutations=2

##################################################################################
# Number of crossover/mating operations to perform each generation.
##################################################################################

num_crossovers=2

##################################################################################
# Size of the population.
##################################################################################

pop_size=5
