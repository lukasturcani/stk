import warnings, os, shutil, sys
warnings.filterwarnings("ignore")

from .ga import (Population, GATools,
                 GAInput, InputHelp, Normalization)
from .convenience_tools import (time_it, tar_output,
                                archive_output, kill_macromodel)
from .ga import plotting as plot

def print_info(info):
    """
    Prints `info` and underlines it.

    """
    print('\n\n' + info + '\n' + '-'*len(info), end='\n\n')

def run():
    """
    Runs the GA.

    """

    # Save the current directory as the `launch_dir`.
    launch_dir = os.getcwd()

    # Running MacroModel optimizations sometimes leaves applications
    # open.This closes them. If this is not done, directories may not
    # be possible to move.
    kill_macromodel()

    # If an output folder of MMEA exists, archive it. This just moves
    # any ``output`` folder in the cwd to the ``old_output`` folder.
    archive_output()

    # Create a new output directory and move into it. Save its path as
    # the root directory.
    os.mkdir('output')

    # Copy the input script into the output folder - this is useful for
    # keeping track of what input was used to generate the output.
    shutil.copyfile(sys.argv[1], os.path.join('output',
                                     os.path.split(sys.argv[1])[-1]))

    os.chdir('output')
    root_dir = os.getcwd()
    os.mkdir('counters')
    # Make template names for the counters which are created showing
    # which members get selected.
    mutation_counter = os.path.join(root_dir, 'counters',
                            "gen_{}_mutation_counter.png")
    crossover_counter = os.path.join(root_dir, 'counters',
                                "gen_{}_crossover_counter.png")
    gen_counter = os.path.join(root_dir, 'counters',
                                "gen_{}_selection_counter.png")

    # Get the name of the input file and load its contents into a
    # ``GAInput`` instance. Info about input file structure is
    # documented in ``GAInput`` docstring.
    ga_input = GAInput(os.path.basename(sys.argv[1]))
    # Load all molecules stored in databases into memory.
    for db in ga_input.databases:
        Population.load(db, load_names=False)

    # Make a Population which stores all previous generations to keep
    # track of progress.
    progress = Population(ga_input.ga_tools())
    # Make a Population which stores every molecule, selected or not,
    # made during the GA run.
    run_db = Population()
    # The variable `id_` is used to give each molecule a unique name
    # during the GA run.
    id_ = 0

    # Make a ``pop_dumps`` directory for string populations as the GA
    # progresses. For debugging, allows the restoration of state before
    # a crash.
    os.mkdir('pop_dumps')
    pop_dump_path = os.path.join(root_dir, 'pop_dumps')
    # Make the ``scratch`` directory which acts as the working
    # directory during the GA run.
    os.mkdir('scratch')
    os.chdir('scratch')

    # Generate the initial population.
    with time_it():
        pop_init = getattr(Population, ga_input.initer().name)
        print_info('Generating initial population.')

        # If the initialization function is ``load()`` a restart run is
        # assumed.
        if pop_init.__name__ == 'load':
            pop = pop_init(**ga_input.initer().params,
                           ga_tools=ga_input.ga_tools())
            ga_input.pop_size = len(pop)

            for mem in pop:
                name = os.path.basename(mem.file)

                mem.file = os.path.join(os.getcwd(), name)

        else:
            pop = pop_init(**ga_input.initer().params,
                           size=ga_input.pop_size,
                           ga_tools=ga_input.ga_tools())

    # Give each population member a name for easy identification when
    # logging.
    for mem in pop:
        mem.name = str(id_)
        id_ += 1

    # Dump the population before attempting any operations.
    pop.dump(os.path.join(pop_dump_path, 'init_pop.json'))

    with time_it():
        print_info('Optimizing the population.')
        pop = Population(pop.ga_tools, *pop.optimize_population())

    with time_it():
        print_info('Calculating the fitness of population members.')
        pop = Population(pop.ga_tools, *pop.calculate_member_fitness())

    with time_it():
        print_info('Normalizing fitness values.')
        pop.normalize_fitness_values()

    # Print the scaled and unscaled fitness values.
    for macro_mol in sorted(pop, reverse=True):
        print(macro_mol.name)
        print(macro_mol.fitness, '-', macro_mol.unscaled_fitness)
        print('\n')

    # Save the generation.
    with time_it():
        print_info('Recording progress.')
        progress.add_subpopulation(pop)
        run_db.add_members(pop)

    # Run the GA.
    for x in range(1, ga_input.num_generations+1):
        # Check that the population has the correct size.
        assert len(pop) == ga_input.pop_size

        print_info('Generation {} of {}.'.format(x,
                                             ga_input.num_generations))

        with time_it():
            print_info('Starting crossovers.')
            offspring = pop.gen_offspring(crossover_counter.format(x))

        with time_it():
            print_info('Starting mutations.')
            mutants = pop.gen_mutants(mutation_counter.format(x))

        with time_it():
            print_info('Adding offsping and mutants to population.')
            pop += offspring + mutants

        with time_it():
            print_info('Removing duplicates, if any.')
            pop.remove_duplicates()

        # Make sure that every population member has a name.
        for mem in pop:
            if not mem.name:
                mem.name = str(id_)
                id_ += 1

        # Dump the population before attempting any operations.
        pop.dump(os.path.join(pop_dump_path,
                              'gen_{}_unselected.json'.format(x)))

        with time_it():
            print_info('Optimizing the population.')
            pop = Population(pop.ga_tools, *pop.optimize_population())

        with time_it():
            print_info(('Calculating the fitness'
                        ' of population members.'))
            pop = Population(pop.ga_tools,
                             *pop.calculate_member_fitness())

        with time_it():
            print_info('Normalizing fitness values.')
            pop.normalize_fitness_values()

        # Print the scaled and unscaled fitness values.
        for macro_mol in sorted(pop, reverse=True):
            print(macro_mol.name)
            print(macro_mol.fitness, '-', macro_mol.unscaled_fitness)
            print('\n')

        with time_it():
            print_info('Selecting members of the next generation.')
            pop = pop.gen_next_gen(ga_input.pop_size,
                                   gen_counter.format(x))

        # Save the generation.
        with time_it():
            print_info('Recording progress.')
            progress.add_subpopulation(pop)
            run_db.add_members(pop)
            pop.dump(os.path.join(pop_dump_path,
                                  'gen_{}_selected.json'.format(x)))

        # If the user defined some premature exit function, check if
        # the exit criterion has been fulfilled.
        if pop.exit():
            break

    # Running MacroModel optimizations sometimes leaves applications
    # open. This closes them. If this is not done, directories may not
    # be possible to move.
    kill_macromodel()

    os.chdir(root_dir)
    # Dump the `progress` and `run_db` populations.
    progress.dump('progress.json')
    run_db.dump('run_db.json')

    # Plot the results of the GA run.
    with time_it():
        print_info('Plotting EPP.')
        # Remove any failed molecules.
        progress.remove_failures()
        # Make sure all fitness values are normalized.
        progress.normalize_fitness_values()

        plot.fitness_epp(progress, 'epp.png')
        plot.parameter_epp(progress, 'epp.png')

    # Remove the ``scratch`` directory.
    shutil.rmtree('scratch')
    # Write the .mol files of the final population.
    pop.write('final_pop', True)

    # Move the ``output`` folder into the ``old_output`` folder.
    os.chdir(launch_dir)

    with time_it():
        print_info('Compressing output.')
        tar_output()

    with time_it():
        archive_output()

def helper():
    """
    Takes care of the -h option.

    """

    InputHelp(sys.argv[-1])

def compare():
    """
    Takes care of the -c option.

    """

    launch_dir = os.getcwd()

    # If an output folder of MMEA exists, archive it. This moves any
    # ``output`` folder in the cwd to the ``old_output`` folder.
    archive_output()

    # Create a new output directory.
    os.mkdir('output')

    # Copy the input script into the output folder - this is useful for
    # keeping track of what input was used to generate the output.
    shutil.copyfile(sys.argv[2], os.path.join('output',
                                   os.path.split(sys.argv[2])[-1]))

    # Get the fitness and normaliztion function data from the input
    # file.
    inp = GAInput(sys.argv[2])

    # Create the encapsulating population.
    pop = Population()
    # Load the fitness and normalization functions into the population.
    pop.ga_tools.input = inp
    pop.ga_tools.fitness = inp.fitness_func
    pop.ga_tools.normalization = (
                            Normalization(inp.normalization_func) if
                            inp.normalization_func else None)

    # Load the populations you want to compare, calculate the fitness
    # of members and place them into the encapsulating population.
    os.chdir('output')
    for i, pop_path in enumerate(inp.comparison_pops):
        sp_dir = os.path.join(os.getcwd(), 'pop{}'.format(i))
        sp = Population.load(pop_path, GATools.init_empty())
        sp.ga_tools.fitness = pop.ga_tools.fitness
        sp.ga_tools.normalization = pop.ga_tools.normalization

        # Before calculating fitness, remove data from previous fitness
        # calculations. Also update the file name with the new
        # directory.
        for ind in sp:
            _, name = os.path.split(ind.file)
            ind.file = os.path.join(sp_dir, name)
            ind.unscaled_fitness = None
            ind.fitness_fail = True
            ind.fitness = None
            ind.progress_params = None

        sp.write(sp_dir)

        # Calculation of fitness is done here, not in the overall pop,
        # to maintain structure.
        sp = Population(*sp.calculate_member_fitness(), sp.ga_tools)
        pop.add_subpopulation(sp)

    # Only try to run a normalization function if one was defined in
    # the input file.
    if inp.normalization_func:
        pop.normalize_fitness_values()
        plot.progress_params(pop, 'param_comparison.png')

    plot.subpopulations(pop, 'fitness_comparison.png')

    # Print the scaled and unscaled fitness values.
    for macro_mol in sorted(pop, reverse=True):
        print(macro_mol.file)
        print(macro_mol.fitness, '-', macro_mol.unscaled_fitness)
        print('\n')

    os.chdir(launch_dir)
    archive_output()

if __name__ == '__main__':
    if '-h' in sys.argv:
        helper()
    elif '-c' in sys.argv:
        compare()
    else:
        run()
