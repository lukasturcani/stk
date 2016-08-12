import os
import shutil
import sys
from itertools import islice

from .classes import (Population, GATools, Selection, Mutation, Mating, 
                      FunctionData, FourPlusSix, EightPlusTwelve, 
                      GAInput)
                      
# Get the name of the input file and load its contents into a 
# ``GAInput`` instance. Info about input file structure is documented in
# ``GAInput`` docstring.
ga_input = GAInput(sys.argv[1])

# If an output folder of MMEA exists, archive it.
if 'output' in os.listdir():
    if 'old_output' not in os.listdir():
        os.mkdir('old_output')
    num = len(os.listdir('old_output'))
    new_dir = os.path.join('old_output', str(num))
    print('Moving old output dir.')
    shutil.copytree('output', new_dir)
    shutil.rmtree('output')
    
# Create a new output directory and move into it. Save its path as the
# root directory.
os.mkdir('output')
os.chdir('output')
root_dir = os.getcwd()
os.mkdir('initial')
os.chdir('initial')

# Use data from the input file to create a ``GATools`` instance for the
# main population.
selector = Selection(ga_input.generational_select_func, 
                     ga_input.parent_select_func, 
                     ga_input.mutant_select_func)
mator = Mating(ga_input.mating_func, ga_input.num_matings)
mutator = Mutation(ga_input.mutation_func, ga_input.num_mutations)
ga_tools = GATools(selector, mator, mutator, ga_input.opt_func)

# Generate an initial population.
pop_init = getattr(Population, ga_input.init_func.name)
pop = pop_init(**ga_input.init_func.params, 
               size=ga_input.pop_size, 
               ga_tools=ga_tools)
# Run the GA.
for x in range(ga_input.num_generations):
    print(('Generation {0} started. Stop at generation {1}. '
            'Population size is {2}.').format(x, 
                                             ga_input.num_generations-1, 
                                             len(pop)))
    # At the start of each generation go into the root directory and 
    # create a folder to hold the next generation's ``.mol`` files.
    # Change into the newly created directory.
    os.chdir(root_dir)
    os.mkdir(str(x))
    os.chdir(str(x))
    
    print('Staring mating.')
    offspring = pop.gen_offspring()
    
    print('Starting mutations.')
    mutants = pop.gen_mutants()
    
    print('Adding offsping and mutants to population.')
    pop += offspring + mutants
    
    print('Optimizing population')
    pop.optimize_population()
    
    print('Selecting members of the next generation.')
    pop = Population(ga_tools, *(islice(pop.select('generational'),
                                        0, ga_input.pop_size)))

    # Create a folder within a generational folder for the the ``.mol``
    # files corresponding to molecules selected for the next generation.
    # Place the ``.mol`` files into that folder.
    print('Copying .mol files to population directory.')
    os.mkdir('selected')
    os.chdir('selected')
    pop.write_population_to_dir(os.getcwd())
    
        