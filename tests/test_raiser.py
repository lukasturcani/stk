# Import stuff.
import mtk
import os
from os.path import join
from glob import iglob

# ####################################################################
# Define variables.
# ####################################################################

cx1 = '/work/lt912/mtk/dbs/liverpool_refined'
local = '/Users/Enrico/Documents/PostDoc/Database/liverpool_refined/'

path = next(x for x in [cx1, local] if os.path.exists(x))

bb_db = join(path, 'aldehydes_3f')
lk_db = join(path, 'amines_2f')


databases = iglob(join(path, '*.json'))

# ####################################################################
# Population initialization function.
# ####################################################################

lks = [mtk.StructUnit2(filename) for filename in iglob(join(lk_db, '*.mol'))]
bbs = [mtk.StructUnit3(filename) for filename in iglob(join(bb_db, '*.mol'))]

init_func = {'NAME': 'init_random',
             'macromol_class': mtk.Cage,
             'building_blocks': [bbs, lks],
             'topologies': [mtk.FourPlusSix()]}

pop_size = 10

# ####################################################################
# Selection function for selecting the next generation.
# ####################################################################

generational_select_func = {'NAME': 'stochastic_sampling',
                            'use_rank': True}

# ####################################################################
# Selection function for selecting parents.
# ####################################################################

crossover_select_func = {'NAME': 'crossover_roulette'}

# ####################################################################
# Selection function for selecting molecules for mutation.
# ####################################################################

mutation_select_func = {'NAME': 'stochastic_sampling',
                        'duplicates': True}

# ####################################################################
# Crossover functions.
# ####################################################################

crossover_funcs = [{'NAME': 'bb_lk_exchange'}]

# ####################################################################
# Mutation functions.
# ####################################################################

mutation_func1 = {'NAME': 'random_bb',
                  'mols': bbs,
                  'key': lambda x: x.__class__ is mtk.StructUnit3}

mutation_func2 = {'NAME': 'random_bb',
                  'mols': lks,
                  'key': lambda x: x.__class__ is mtk.StructUnit2}

mutation_func3 = {'NAME': 'similar_bb',
                  'mols': bbs,
                  'key': lambda x: x.__class__ is mtk.StructUnit3}

mutation_func4 = {'NAME': 'similar_bb',
                  'mols': lks,
                  'key': lambda x: x.__class__ is mtk.StructUnit2}

mutation_funcs = [mutation_func1,
                  mutation_func2,
                  mutation_func3,
                  mutation_func4]

# ####################################################################
# Optimization function.
# ####################################################################

opt_func = {'NAME': 'raiser',
            'param1': 1,
            'param2': 2}

# ####################################################################
# Fitness function.
# ####################################################################

fitness_func = {'NAME': 'raiser',
                'param1': 1,
                'param2': 2}

# ####################################################################
# Normalization functions.
# ####################################################################

# First shift all energy values so that they are always positive.
norm_func1 = {'NAME': 'shift_elements',
              'indices': [-1]}

# Second, make sure that the magnitudes of all fitness parameters
# are comparable.
norm_func2 = {'NAME': 'magnitudes'}

# Third, combine the fitness parameter into a single fitness value.
norm_func3 = {'NAME': 'combine',
              'coefficients': [1, 1, 1, 1],
              'exponents': [1, 1, 1, 1]}

# Last, invert the fitness value because all fitness parameters
# are inversely proportional to fitness.
norm_func4 = {'NAME': 'invert'}

# Now define the normalization_funcs variable holding all of the
# normalization functions to be used.

normalization_funcs = [norm_func1, norm_func2, norm_func3, norm_func4]

# ####################################################################
# Number of generations to create.
# ####################################################################

num_generations = 15

# ####################################################################
# Number of mutation operations to perform each generation.
# ####################################################################

num_mutations = 3

# ####################################################################
# Number of crossover/mating operations to perform each generation.
# ####################################################################

num_crossovers = 3
