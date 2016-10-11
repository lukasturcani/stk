import os
import shutil
import sys
import warnings

from .classes import (Population, GATools, Selection, Mutation, Mating, 
                      GAInput)
from .classes.exception import PopulationSizeError
from .optimization import kill_macromodel
from .convenience_functions import time_it

# Running MacroModel optimizations sometimes leaves applications open.
# This closes them. If this is not done, directories may not be possible
# to move. 
kill_macromodel()
                      
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

    # Wait for the copy to complete before removing the old folder.
    mv_complete = False    
    while not mv_complete:
        try:
            shutil.rmtree('output')
            mv_complete = True
        except:
            pass
    
# Create a new output directory and move into it. Save its path as the
# root directory.
    
# Wait for previous operations to finish before starting these.
mk_complete = False    
while not mk_complete:
    try:
        os.mkdir('output')
        mk_complete = True
    except:
        pass    
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
mutator = Mutation(ga_input.mutation_func, ga_input.num_mutations,
                   weights=ga_input.mutation_weights)
ga_tools = GATools(selector, mator, mutator, 
                   ga_input.opt_func, ga_input.fitness_func)
ga_tools.ga_input = ga_input

# Generate and optimize an initial population.
with time_it():
    pop_init = getattr(Population, ga_input.init_func.name)
    print(('\n\nGenerating initial population.\n'
         '------------------------------\n\n'))
    pop = pop_init(**ga_input.init_func.params, 
                   size=ga_input.pop_size, 
                   ga_tools=ga_tools)

with time_it():    
    print(('\n\nOptimizing the population.\n'
          '--------------------------\n\n'))
    pop = Population(pop.ga_tools, *pop.optimize_population())

with time_it():    
    print('\n\nCalculating the fitness of population members.\n'
        '----------------------------------------------\n\n') 
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pop.calculate_member_fitness()

# Run the GA.
for x in range(ga_input.num_generations):
    # Save the min, max and mean values of the population.    
    pop.progress_update()
    
    # Check that the population has the correct size.
    if len(pop) != ga_input.pop_size:
        raise PopulationSizeError('Population has the wrong size.')
    
    print(('\n\nGeneration {0} started. Stop at generation {1}. '
            'Population size is {2}.\n'
            '---------------------------------------------------------'
            '--------------\n\n').format(x, ga_input.num_generations-1, 
                                            len(pop)))
    # At the start of each generation go into the root directory and 
    # create a folder to hold the next generation's ``.mol`` files.
    # Change into the newly created directory.
    os.chdir(root_dir)
    os.mkdir(str(x))
    os.chdir(str(x))
    
    with time_it():
        print('\n\nStarting mating.\n----------------\n\n')
        offspring = pop.gen_offspring()

    with time_it():
        print('\n\nStarting mutations.\n-------------------\n\n')
        mutants = pop.gen_mutants()
    
    with time_it():
        print(('\n\nAdding offsping and mutants to population.'
              '\n------------------------------------------\n\n'))
        pop += offspring + mutants
    
    with time_it():
        print(('\n\nRemoving duplicates, if any.\n'
               '----------------------------\n\n')    )
        pop.remove_duplicates()        
    
    pop.dump(os.path.join(os.getcwd(), 'pop_dump'))    
    
    with time_it():        
        print(('\n\nOptimizing the population.\n'
              '--------------------------\n\n'))
        pop = Population(pop.ga_tools, *pop.optimize_population())

    with time_it():        
        print('\n\nCalculating the fitness of population members.\n'
            '----------------------------------------------\n\n')    
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pop.calculate_member_fitness()    

    with time_it():        
        print(('\n\nSelecting members of the next generation.\n'
               '-----------------------------------------\n\n'))
        pop = pop.gen_next_gen(ga_input.pop_size)
        
    # Create a folder within a generational folder for the the ``.mol``
    # files corresponding to molecules selected for the next generation.
    # Place the ``.mol`` files into that folder.
    print(('\n\nCopying .mol files to population directory.\n'
           '-------------------------------------------\n\n'))
    os.mkdir('selected')
    os.chdir('selected')
    pop.write_population_to_dir(os.getcwd())
    
# Running MacroModel optimizations sometimes leaves applications open.
# This closes them. If this is not done, directories may not be possible
# to move.     
kill_macromodel()

# Update a final time and plot the results of the GA run.
pop.progress_update()
pop.plot_epp(os.path.join(root_dir, 'epp.png'))

