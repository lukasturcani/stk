#############################################################
# Define basic GA constants.
#############################################################

num_generations = 250

# Number of mutation operations to perform each generation.
num_mutations = 10

# Number of crossover operations to perform each generation.
num_crossovers = 5

# Size of the population.
pop_size = 25

#############################################################
# Define variables.
#############################################################

monomer_database = '/path/to/monomer/database'
polymer_db1 = 'path/to/previously/made/polymers1.json'
polymer_db2 = '/path/to/previously/made/polymers2.json'

#############################################################
# Databases of stored molecules to use.
#############################################################

databases = [polymer_db1, polymer_db1]

#############################################################
# Population initialization function.
#############################################################

init_func =  {

              'NAME' : 'init_random_polymers',
              'monomer_db' : monomer_database,
              'topologies' : [Linear('AB', [0, 1], 24),
                              Linear('ABBC', [0, 1, -1, 0], 12)]
            }

#############################################################
# Selection function for selecting the next generation.
#############################################################

generational_select_func =  {
                             'NAME' : 'stochastic_sampling',
                             'use_rank' : True
                            }

#############################################################
# Selection function for selecting parents.
#############################################################

parent_select_func =  {'NAME' :'crossover_roulette'}

#############################################################
# Selection function for selecting molecules for mutation.
#############################################################

mutant_select_func =  {
                       'NAME' : 'stochastic_sampling',
                       'duplicates':True
                      }

#############################################################
# Crossover functions.
#############################################################

# Crossover function 1.
crossover_func1 =  {'NAME' : 'polymer_monomer_shuffle'}

# Crossover function 2.
crossover_func2 =  {'NAME' : 'polymer_topology_shuffle'}

# Tell MMEA which crossover functions to use.
crossover_funcs = [crossover_func1, crossover_func2]

# Probability that each crossover function is selected.
crossover_weights = [0.25, 0.75]

#############################################################
# Mutation functions
#############################################################

# Mutation function 1.
mutation_func1 =  {
                   'NAME' : 'cage_random_bb',
                   'database' : bb_db_path
                  }

# Mutation function 2.
mutation_func2 = {
                  'NAME' : 'cage_random_lk',
                  'database' : lk_db_path
                 }

# Tell MMEA which mutation functions to use.
mutation_funcs = [mutation_func1, mutation_func2]

# Probability that each mutation function is selected.
mutation_weights = [0.5, 0.5]

#############################################################
# Optimization function.
#############################################################

opt_func =  {
              'NAME' : 'rdkit_optimization',
              'embed' : True
            }

#############################################################
# Fitness function.
#############################################################

# Calculates the polymer's energy and surface area.
fitness_func =  {'NAME' : 'polymer_E_and_SA'}

#############################################################
# Normalization functions.
#############################################################

# First shift all energy values so that they are always positive.
norm_func1 = {
              'NAME' : 'shift_elements',
              'indices' : [0]
             }

# Second, make sure that the magnitudes of all fitness parameters
# are comparable.
norm_func2 =  {'NAME' : 'magnitudes'}

# Combine the fitness parameters into a single fitness value.
norm_func3 = {
              'NAME' : 'combine',
              'coefficients' : [1, 3],
              'exponents' : [1, 1]
             }

normalization_funcs = [norm_func1, norm_func2, norm_func3]
