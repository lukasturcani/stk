import warnings, os, shutil, sys, logging
from os.path import join, basname, abspath
warnings.filterwarnings("ignore")

from .ga import (Population, GATools,
                 GAInput, InputHelp, Normalization)
from .convenience_tools import (time_it, tar_output,
                                archive_output, kill_macromodel)
from .ga import plotting as plot


class ProgressHandler(logging.Handler):
    ...

class DBHandler(logging.Handler):
    ...


logger = logging.getLogger(__name__)


def run():
    """
    Runs the GA.

    """

    # 1. Set up the directory structure.

    *_, ifile = sys.argv
    ifile = abspath(ifile)
    launch_dir = os.getcwd()
    # Running MacroModel optimizations sometimes leaves applications
    # open.This closes them. If this is not done, directories may not
    # be possible to move.
    kill_macromodel()
    # Move any ``output`` dir in the cwd into ``old_output``.
    archive_output()
    os.mkdir('output')
    os.chdir('output')
    root_dir = os.getcwd()
    # If logging level is DEBUG or less, make population dumps as the
    # GA progresses and save images of selection counters.
    mcounter = ccounter = gcounter = None
    if logger.isEnabledFor(logging.DEBUG):
        os.mkdir('counters')
        os.mkdir('pop_dumps')
        cstr = join(root_dir, 'counters', 'gen_{}_')
        mcounter = cstr + 'mutation_counter.png'
        ccounter = cstr + 'crossover_counter.png'
        gcounter = cstr + 'selection_counter.png'
    # Copy the input script into the ``output`` folder.
    shutil.copyfile(ifile, basename(ifile))
    # Make the ``scratch`` directory which acts as the working
    # directory during the GA run.
    os.mkdir('scratch')
    os.chdir('scratch')

    # 2. Initialize the population.

    ga_input = GAInput(ifile)
    logger.info('Loading molecules from any provided databases.')
    for db in ga_input.databases:
        Population.load(db)

    logger.info('Generating initial population.')
    init_func = getattr(Population, ga_input.initer().name)
    pop = init_func(**ga_input.initer().params,
                      size=ga_input.pop_size,
                      ga_tools=ga_input.ga_tools())
    # The variable `id_` is used to give each molecule a unique name
    # during the GA run. Assuming they are not loaded from a database
    # with a name.
    id_ = 0
    for mem in pop:
        if not mem.name:
            mem.name = str(id_)
            id_ += 1

    # Dump the population before attempting any operations.
    pop.dump(join(pop_dump_path, 'init_pop.json'))
    logger.info('Optimizing the population.')
    pop.optimize_population()
    print_info('Calculating the fitness of population members.')
    pop.calculate_member_fitness()
    logger.info('Normalizing fitness values.')
    pop.normalize_fitness_values()

    # Print the scaled and unscaled fitness values.
    for macro_mol in sorted(pop, reverse=True):
        print(macro_mol.name)
        print(macro_mol.fitness, '-', macro_mol.unscaled_fitness)
        print('\n')

    logger.info('Recording progress.')
    progress.add_subpopulation(pop)
    run_db.add_members(pop)

    # 3. Run the GA.

    for x in range(1, ga_input.num_generations+1):
        # Check that the population has the correct size.
        assert len(pop) == ga_input.pop_size
        logger.info('Generation {} of {}.'.format(x, ga_input.num_generations))
        logger.info('Starting crossovers.')
        offspring = pop.gen_offspring(crossover_counter.format(x))
        logger.info('Starting mutations.')
        mutants = pop.gen_mutants(mutation_counter.format(x))
        logger.info('Adding offsping and mutants to population.')
        pop += offspring + mutants
        logger.info('Removing duplicates, if any.')
        pop.remove_duplicates()

        # Make sure that every population member has a name.
        for mem in pop:
            if not mem.name:
                mem.name = str(id_)
                id_ += 1

        # Dump the population before attempting any operations.
        pop.dump(join(pop_dump_path, 'gen_{}_unselected.json'.format(x)))

        logger.info('Optimizing the population.')
        pop.optimize_population()
        logger.info('Calculating the fitness of population members.')
        pop.calculate_member_fitness()
        logger.info('Normalizing fitness values.')
        pop.normalize_fitness_values()

        # Print the scaled and unscaled fitness values.
        for macro_mol in sorted(pop, reverse=True):
            print(macro_mol.name)
            print(macro_mol.fitness, '-', macro_mol.unscaled_fitness)
            print('\n')

        logger.info('Selecting members of the next generation.')
        pop = pop.gen_next_gen(ga_input.pop_size, gen_counter.format(x))
        logger.info('Recording progress.')
        progress.add_subpopulation(pop)
        run_db.add_members(pop)
        pop.dump(join(pop_dump_path, 'gen_{}_selected.json'.format(x)))

        # If the user defined some premature exit function, check if
        # the exit criterion has been fulfilled.
        if pop.exit():
            break

    kill_macromodel()
    os.chdir(root_dir)
    # Dump the `progress` and `run_db` populations.
    progress.dump('progress.json')
    run_db.dump('run_db.json')

    # Plot the results of the GA run.
    logger.info('Plotting EPP.')
    # Make sure all fitness values are normalized.
    progress.normalize_fitness_values()
    plot.fitness_epp(progress, 'epp.png')
    progress.remove_members(lambda x :
          progress.ga_tools.fitness.name not in x.progress_params)
    plot.parameter_epp(progress, 'epp.png')

    # Remove the ``scratch`` directory.
    shutil.rmtree('scratch')
    # Write the .mol files of the final population.
    pop.write('final_pop', True)
    os.chdir(launch_dir)
    logger.info('Compressing output.')
    tar_output()
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
    shutil.copyfile(sys.argv[2], join('output',
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
        sp_dir = join(os.getcwd(), 'pop{}'.format(i))
        sp = Population.load(pop_path, GATools.init_empty())
        sp.ga_tools.fitness = pop.ga_tools.fitness
        sp.ga_tools.normalization = pop.ga_tools.normalization

        # Before calculating fitness, remove data from previous fitness
        # calculations. Also update the file name with the new
        # directory.
        for ind in sp:
            _, name = os.path.split(ind.file)
            ind.file = join(sp_dir, name)
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
