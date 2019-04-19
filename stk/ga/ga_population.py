"""
Defines :class:`GAPopulation`.

:class:`GAPopulation` extends :class:`.Population` by addings tools for
carrying out GA operations.

"""

from collections import Counter
import psutil

from .plotting import plot_counter
from ..population import Population


class GAPopulation(Population):
    """
    Used for applying GA operations to molecules.

    The GA is invoked by calling a number of methods of this class,
    such as :meth:`gen_offspring`, :meth:`gen_mutants` and
    :meth:`select`. However, this class only implements container
    related functionality. It delegates GA operations to the
    :class:`.Crossover`, :class:`.Mutation` and :class:`.Selection`
    classes.

    These classes are organised as follows. Each :class:`GAPopulation`
    instance has a :attr:`ga_tools` attribute. This holds a
    :class:`.GATools` instance. The :class:`.GATools` instance is just
    a container. It holds a :class:`.Crossover`, :class:`.Mutation`
    and :class:`.Selection` instance. These instances deal with the
    :class:`GAPopulation` instance they are held by and perform the
    various crossover, mutation and selection operations on it.
    Any functionality related to the GA should be delegated to these
    instances. The :meth:`gen_offspring` and :meth:`gen_mutants`
    methods can serve as a guide to how this should be done.

    """

    def set_ga_tools(self,
                     generation_selector,
                     mutation_selector,
                     crossover_selector,
                     mutator,
                     crosser,
                     fitness_calculator,
                     fitness_normalizer):
        """

        """

        self.generation_selector = generation_selector
        self.mutation_selector = mutation_selector
        self.crossover_selector = crossover_selector
        self.mutator = mutator
        self.crosser = crosser
        self.fitness_calculator = fitness_calculator
        self.fitness_normalizer = fitness_normalizer

    def calculate_member_fitness(self, processes=psutil.cpu_count()):
        """
        Applies the fitness function to all members.

        The fitness function is defined in the attribute
        :attr:`.GATools.fitness` of the :class:`.GATools` instance
        held in the :attr:`ga_tools` attribute of the population.

        The calculation will be performed serially or in parallel
        depending on the flag :attr:`.GATools.parallel`. The serial
        version may be faster in cases where all molecules have already
        had their fitness values calcluated. This is because all
        calculations will be skipped. In this case creating a parallel
        process pool creates unncessary overhead.

        Notes
        -----
        This method will result in the modification of the
        :attr:`.MacroMolecule.unscaled_fitness` attribute of molecules
        held by the population.

        Parameters
        ----------
        processes : :class:`int`
            The number of parallel processes to create.

        Returns
        -------
        None : :class:`NoneType`

        """

        if processes == 1:
            _calc_fitness_serial(self.ga_tools.fitness, self)
        else:
            _calc_fitness(self.ga_tools.fitness, self, processes)

    def gen_mutants(self):
        """
        Returns a :class:`GAPopulation` of mutants.

        The selection function which decides which molecules are
        selected for mutation is defined in the
        :attr:`.Selection.mutation` attribute of the
        :class:`.Selection` instance held by the population via its
        :attr:`ga_tools` attribute.

        The mutation function(s) to be used are defined in the
        :class:`.Mutation` instance held by the population via its
        :attr:`ga_tools` attribute.

        Parameters
        ----------
        counter_name : :class:`str`, optional
            The name of the ``.png`` file which holds a graph showing
            which members were selected for mutation.

        Returns
        -------
        :class:`GAPopulation`
            A population holding mutants created by mutating the
            :class:`.MacroMolecule` instances held by the population.

        """
        """
        Carries out mutation operations on `population`.

        This function selects members of `population` to be mutated
        and mutates them. This goes on until either all possible
        molecules have been mutated or the required number of
        mutation operations have been performed.

        The mutants generated are returned together in a
        :class:`.GAPopulation` object. Any mutants already present in
        `population` are removed.

        Parameters
        ----------
        population : :class:`.GAPopulation`
            The population who's members are to be mutated.

        counter_path : :class:`str`, optional
            The path to the ``.png`` file showing which members were
            selected for mutation. If ``''`` then no file is made.

        Returns
        -------
        :class:`.GAPopulation`
            A population with all the mutants generated held in
            :attr:`~.Population.members`. This does not include mutants
            which correspond to molecules already present in
            `population`.

        """

        counter_name='mutation_counter.png'
        mutant_pop = GAPopulation(ga_tools=population.ga_tools)
        counter = Counter()

        parent_pool = islice(population.select('mutation'),
                             self.num_mutations)
        for i, parent in enumerate(parent_pool, 1):
            logger.info('Mutation number {}. Finish when {}.'.format(
                                          i, self.num_mutations))
            counter.update([parent])
            func_data = np.random.choice(self.funcs, p=self.weights)
            func = getattr(self, func_data.name)
            logger.info(f'Using {func.__name__}.')

            try:
                mutant = func(parent, **func_data.params)

                # If the mutant was retrieved from the cache, log the
                # name.
                if mutant.name:
                    logger.debug(('Mutant "{}" retrieved from '
                                  'cache.').format(mutant.name))

                mutant_pop.members.append(mutant)

            except Exception as ex:
                errormsg = ('Mutation function "{}()" '
                            'failed on molecule "{}".').format(
                             func_data.name, parent.name)
                logger.error(errormsg, exc_info=True)

        mutant_pop -= population

        if counter_path:
            # Update counter with unselected members.
            for member in population:
                if member not in counter.keys():
                    counter.update({member: 0})
            plot_counter(counter, counter_path)

        return mutant_pop

    def gen_offspring(self):
        """
        Returns a :class:`GAPopulation` of offspring molecules.

        The selection function which decides which molecules are
        selected for crossover is defined in the
        :attr:`.Selection.crossover` attribute of the
        :class:`.Selection` instance held by the population via its
        :attr:`ga_tools` attribute.

        The crossover function(s) to be used are defined in the
        :class:`.Crossover` instance held by the population via its
        :attr:`ga_tools` attribute.

        Parameters
        ----------
        counter_name : :class:`str`, optional
            The name of the ``.png`` file showing which members were
            selected for crossover.

        Returns
        -------
        :class:`GAPopulation`
            A population of offspring, created by crossing
            :class:`.MacroMolecule` instances contained in the
            population.

        """
        counter_name='crossover_counter.png'
        """
        Carries out crossover operations on `population`.

        This function selects members of `population` and crosses
        them until either all possible parents have been crossed or the
        required number of successful crossover operations has been
        performed.

        The offspring generated are returned together in a
        :class:`.GAPopulation` instance. Any molecules that are created
        via crossover and match a molecule present in the original
        population are removed.

        Parameters
        ----------
        population : :class:`.GAPopulation`
            The population instance who's members are to be crossed.

        counter_path : :class:`str`, optional
            The name of the ``.png`` file showing which members were
            selected for crossover. If ``''``, then no file is made.

        Returns
        -------
        :class:`.GAPopulation`
            A population with all the offspring generated held in its
            :attr:`~.Population.members` attribute. This does not
            include offspring which correspond to molecules already
            present in `population`.

        """

        offspring_pop = GAPopulation(ga_tools=population.ga_tools)
        counter = Counter()

        parent_pool = islice(population.select('crossover'),
                             self.num_crossovers)
        for i, parents in enumerate(parent_pool, 1):
            logger.info('Crossover number {}. Finish when {}.'.format(
                                           i, self.num_crossovers))
            counter.update(parents)
            # Get the crossover function.
            func_data = np.random.choice(self.funcs, p=self.weights)
            func = getattr(self, func_data.name)
            logger.info(f'Using {func.__name__}.')

            try:
                # Apply the crossover function and supply any
                # additional arguments to it.
                offspring = func(*parents, **func_data.params)

                # Print the names of offspring which have been returned
                # from the cache.
                for o in offspring:
                    if o.name:
                        logger.debug(('Offspring "{}" retrieved '
                                      'from cache.').format(o.name))

                # Add the new offspring to the offspring population.
                offspring_pop.add_members(offspring)

            except Exception:
                errormsg = ('Crossover function "{}()" failed on '
                            'molecules PARENTS.').format(
                            func_data.name)

                pnames = ' and '.join('"{}"'.format(p.name) for
                                      p in parents)
                errormsg = errormsg.replace('PARENTS', pnames)
                logger.error(errormsg, exc_info=True)

        # Make sure that only original molecules are left in the
        # offspring population.
        offspring_pop -= population

        if counter_path:
            # Update counter with unselected members and plot counter.
            for member in population:
                if member not in counter.keys():
                    counter.update({member: 0})
            plot_counter(counter, counter_path)

        return offspring_pop

    def normalize_fitness_values(self):
        """
        Applies the normalization functions on the population.

        Normalization functions scale or modify the fitness values of
        molecules in the population.

        The normalization functions which are applied on the
        population, along with their order, are defined in the
        :class:`.Normalization` instance held by the population via the
        :attr:`ga_tools` attribute.

        Notes
        -----
        This method modifies the :attr:`.MacroMolecule.fitness`
        attribute of molecules held by the population.

        Returns
        -------
        None : :class:`NoneType`

        """

        return
