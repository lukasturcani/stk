"""
Defines selection functions via :class:`Selection`.

.. _`adding selection functions`:

Extending stk: Adding selection functions.
------------------------------------------

If a new selection operation is to be added to ``stk`` it should be
added as a method in :class:`Selection`. The only requirements are
that the first argument is `population` (excluding any `self` or `cls`
arguments).

The naming requirement exists to help users identify which arguments
are handled automatically by the GA and which they need to define in
the input file. The convention is that if the selection function takes
an argument called  `population` it does not have to be specified in
the input file.

All selection functions should be defined as generators, which yield a
member of the population. In the case of crossover selection functions,
they should yield a tuple of members. Generational selection functions
should yield each molecule at most once.

If the selection function does not fit neatly into a single function
make sure that any helper functions are private, i.e. that their names
start with a leading underscore.

"""

import itertools
import numpy as np


class Selection:
    """
    A class for handling all types of selection in the GA.

    Whenever a population needs to have some of its members selected
    for GA operations such as mutation, selection or generational
    selection, it delegates this task to an instance of this class. The
    population holds the instance in :attr:`.GAPopulation.ga_tools`.

    To illustrate how a :class:`Selection` instance works it is best
    to use an example.

    .. code-block:: python

        # 1. Create FunctionData objects which describe the selection
        #    functions to be used. These will correspond to methods
        #    in the Selection class.

        # For selecting members of the next generation.
        gen_select_fn = FunctionData('fittest')

        # For selecting molecules for crossover.
        crossover_select_fn = FunctionData('crossover_roulette', n=5)

        # For selecting molecules to mutate.
        mutation_select_fn = FunctionData('stochastic_sampling',
                                          elites=2, duplicates=True)



        # 2. Create a Selection instance using the chosen functions.

        sel = Selection(gen_select_fn,
                        crossover_select_fn,
                        mutation_select_fn)

        # 3. Create a population from which "sel" will select members.
        pop = GAPopulation(mol1, mol2, mol3, mol4)

        # 4. Now the selection instance can be used to select members
        #    of "pop".

        # To select members of the next generation, call the selection
        # instance with the population as an argument, and a string
        # to indicate generational selection is desired.
        next_gen = sel(pop, 'generational')

        # Note that next_gen is a generator, which yields the
        # selected molecules.
        next_gen  # <generator object ... >

        for mol in next_gen:
            mol  # A molecule in "pop".

        # To select members for crossover, the same idea as for
        # generational selection applies.
        crossover_mols = sel(pop, 'crossover')

        # Again, "crossover_mols" is a generator. This time it yields
        # tuples of molecules however. This is because 2 or more
        # molecules are required per crossover operation.


        # Finally, to select molecules for mutation.
        mols_to_mutate = sel(pop, 'mutation')

        # Again, "mols_to_mutate" is a generator yielding single
        # molecules.

    By providing the strings ``'generational'``, ``'crossover'`` and
    ``'mutation'``, the :class:`Selection` instance knows to pick the
    correct function from the :class:`.FunctionData` instances provided
    during initialization.

    Note that duplicates are not considered during selection. If a
    molecule is present in the population twice, its chance of
    selection does not double.

    Attributes
    ----------
    generational : :class:`.FunctionData`
        This holds the :class:`.FunctionData` object representing the
        selection function used for selecting the next generation of
        molecules.

    crossover : :class:`.FunctionData`
        Holds the :class:`.FunctionData` object representing the
        selection function used for selecting molecules for crossover.

    mutation : :class:`.FunctionData`
        Holds the :class:`.FunctionData` object representing the
        selection function used for selecting molecules to mutate.

    """

    def __init__(self, generational, crossover, mutation):
        """
        Initializes a :class:`Selection` instance.

        Parameters
        ----------
        generational : :class:`.FunctionData`
            This holds the :class:`.FunctionData` object representing
            the selection function used for selecting the next
            generation of molecules.

        crossover : :class:`.FunctionData`
            Holds the :class:`.FunctionData` object representing the
            selection function used for selecting molecules for
            crossover.

        mutation : :class:`.FunctionData`
            Holds the :class:`.FunctionData` object representing the
            selection function used for selecting molecules to mutate.

        """

        self.generational = generational
        self.crossover = crossover
        self.mutation = mutation

    def __call__(self, population, type_):
        """
        Yields members from population.

        Parameters
        ----------
        population : :class:`.GAPopulation`
            The population from which members should be selected.

        type_ : :class:`str`
            The name of a :class:`.Selection` attribute. The name
            corresponds to the type of selection that is desired:
            ``'generational'``, ``'crossover'`` or ``'mutation'``.

        Yields
        -------
        :class:`.MacroMolecule` or :class:`tuple`
            :class:`.MacroMolecule` instances are yielded unless
            `type_` is ``'crossover'``. In this case a tuple of
            :class:`.MacroMolecule` instances is yielded.

        """

        # Make a population without any duplicates.
        unique_pop = population.__class__(*population)
        unique_pop.remove_duplicates()

        # Get the attribute with the name provided as a string in
        # `type_`. This returns a ``FunctionData`` object holding the
        # name of the method which needs to be implemented and any
        # required parameters and their values.
        func_data = getattr(self, type_)
        # Using the name of the method, get the method object with
        # ``getattr``.
        func = getattr(self, func_data.name)
        # Call the method on the population and provide any additional
        # parameters which may be necessary.
        yield from func(unique_pop, **func_data.params)

    def fittest(self, population):
        """
        Yields members of the population, fittest first.

        Parameters
        ----------
        population : :class:`.GAPopulation`
            The population who's members need to be selected.

        Yields
        ------
        :class:`.MacroMolecule`
            The next fittest :class:`.MacroMolecule` instance held by
            `population`.

        """

        for ind in sorted(population, reverse=True):
            yield ind

    def roulette(self, population, elites=0,
                 truncation=None, duplicates=False):
        """
        Yields individuals using roulette selection.

        In roulette selection the probability an individual is selected
        is given by its fitness. If the total fitness is the sum of all
        fitness values, the chance an individuals is selected is given
        by::

            p = individual fitness / total fitness,

        where ``p`` is the probability of selection [#]_.

        Parameters
        ----------
        population : :class:`.GAPopulation`
            The population from which individuals are to be selected.

        elites : :class:`int`, optional
            The number of the fittest members which are guaranteed to
            be yielded first.

        truncation : :class:`int`, optional
            The number of least fit members which will never be
            yielded.

        duplicates : :class:`bool`, optional
            If ``True`` the same member can be yielded more than
            once.

        Yields
        ------
        :class:`.MacroMolecule`
            The next selected population member.

        References
        ----------
        .. [#] http://tinyurl.com/csc3djm

        """

        yielded = set()
        truncation = truncation if truncation is None else -truncation
        pop = sorted(population, reverse=True)[:truncation]

        for x in range(elites):
            yielded.add(pop[x]) if not duplicates else None
            yield pop[x]

        while True:
            valid_pop = [ind for ind in pop if ind not in yielded]

            if not valid_pop:
                break

            total = sum(ind.fitness for ind in valid_pop)
            weights = [ind.fitness / total for ind in valid_pop]

            selected = np.random.choice(valid_pop, p=weights)
            yielded.add(selected) if not duplicates else None
            yield selected

    def deterministic_sampling(self, population, truncation=None):
        """
        Yields individuals using deterministic sampling.

        This algorithm can only be used for selection of mutants. It
        will return duplicates. The results of :meth:`fittest` are
        equivalent to algorithm without duplicates.

        In determnistic sampling the mean fitness value of the
        population is calculated, ``<f>``. For each individual a
        normalized fitness value is then calculated via::

            fn = f / <f>

        where ``fn`` is the normalized fitness value and ``f`` is the
        original fitness value. An individual will be selected ``n``
        times, where ``n`` is the integer value of ``fn``. After all
        individuals where ``n > 0`` are yielded, all inidividuals are
        yielded again, largest decimal part of ``fn`` first.

        Parameters
        ----------
        population : :class:`.GAPopulation`
            The population from which individuals are to be selected.

        truncation : :class:`int`, optional
            The number of least fit members which will never be
            yielded.

        Yields
        ------
        :class:`.MacroMolecule`
            The next selected population member.

        """

        truncation = truncation if truncation is None else -truncation
        pop = sorted(population, reverse=True)[:truncation]
        mean = np.mean([mem.fitness for mem in pop])
        decimals = []
        for mem in pop:
            q, r = divmod(mem.fitness, mean)
            decimals.append((r, mem))
            for i in range(int(q)):
                yield mem

        for r, mem in sorted(decimals, reverse=True):
            yield mem

    def stochastic_sampling(self, population,
                            elites=0, truncation=None,
                            duplicates=False, use_rank=False):
        """
        Yields individuals via stochastic sampling.

        Each fitness value is used to calculate the normalized fitness
        of an individual::

            fn = f / <f>

        where ``fn`` is the normalized fitness value, ``f`` is the
        original fitness value and ``<f>`` is the mean fitness of the
        population. If ``fn`` is greater than ``1`` then the individual
        is guaranteed to be selected. If duplicates are allowed,
        individuals are guaranteed to be selected ``n`` times, where
        ``n`` is the integer part of ``fn``.

        Once all the guarnteed individuals have been yielded, the
        remaining individuals are yielded via the roulette method. The
        weights in the roulette method are based on the decimal part of
        ``fn``.

        Parameters
        ----------
        population : :class:`.GAPopulation`
            The population from which individuals are to be selected.

        elites : :class:`int`, optional
            The number of the fittest members which are guaranteed to
            be yielded first.

        truncation : :class:`int`, optional
            The number of least fit members which will never be
            yielded.

        duplicates : :class:`bool`, optional
            If ``True`` the same member can be yielded more than
            once.

        use_rank : :class:`bool`, optional
            When ``True`` the fitness value of an individual is
            calculated as, ``f = 1/rank``.

        Yields
        ------
        :class:`.MacroMolecule`
            The next selected population member.

        """

        yielded = set()
        truncation = truncation if truncation is None else -truncation
        pop = sorted(population, reverse=True)[:truncation]

        for x in range(elites):
            yielded.add(pop[x])
            yield pop[x]

        for r, mem in enumerate(pop, 1):
            mem._ssfitness = 1 / r if use_rank else mem.fitness

        yield from self._stochastic_sampling_guaranteed(pop,
                                                        yielded,
                                                        duplicates)
        yield from self._stochastic_sampling_roulette(pop,
                                                      yielded,
                                                      duplicates)

    def _stochastic_sampling_guaranteed(self, pop,
                                        yielded, duplicates):
        """
        Yielded the members guaranteed by stochastic sampling.

        Parameters
        ----------
        pop : :class:`list` of :class:`.MacroMolecule`
            The molecules which are being selected.

        yielded : :class:`set` of :class:`.MacroMolecule`
            Holds all previously yielded molecules.

        duplicates : :class:`bool`
            Indicates whether a member can be yielded more than once.

        Yields
        ------
        :class:`.MacroMolecule`
            The next member whose normalized fitness integer component
            is greater than 0.

        """

        mean = np.mean([mem._ssfitness for mem in pop])
        for mem in pop:
            q = int(divmod(mem._ssfitness, mean)[0])
            # Account for the fact that a molecule may have been
            # yielded due to elitism.
            q -= 1 if mem in yielded else 0
            q = 0 if mem in yielded and not duplicates else q
            for x in range(q):
                yielded.add(mem)
                yield mem
                if not duplicates:
                    break

    def _stochastic_sampling_roulette(self, pop,
                                      yielded, duplicates):
        """
        Does the roulette component of stochastic sampling.

        Parameters
        ----------
        pop : :class:`list` of :class:`.MacroMolecule`
            The molecules which are being selected.

        yielded : :class:`set` of :class:`.MacroMolecule`
            Holds all previously yielded molecules.

        duplicates : :class:`bool`
            Indicates whether a member can be yielded more than once.

        Yields
        ------
        :class:`.MacroMolecule`
            The seleceted population member.

        """

        yielded = set() if duplicates else yielded

        while True:
            valid_pop = [ind for ind in pop if ind not in yielded]

            if not valid_pop:
                break

            mean = np.mean([mem._ssfitness for mem in pop])
            rtotal = sum(divmod(ind._ssfitness, mean)[1] for
                         ind in valid_pop)
            weights = [divmod(ind._ssfitness, mean)[1] / rtotal for
                       ind in valid_pop]

            selected = np.random.choice(valid_pop, p=weights)
            yielded.add(selected) if not duplicates else None
            yield selected

    def all_combinations(self, population):
        """
        Yields every possible pairing of individuals from a population.

        This yields members regardless of which subpopulation they are
        in. Each pair is only returned once. This means if ``(1,2)`` is
        yielded ``(2,1)`` will not be.

        Parameters
        ----------
        population : :class:`.GAPopulation`
            The population from which parents should be selected.

        Yields
        ------
        :class:`tuple` of :class:`.MacroMolecule`
            The :class:`.MacroMolecule` instances which are to be
            crossed.

        """

        # Get all combinations of size 2.
        for mol1, mol2 in itertools.combinations(population, 2):
            yield mol1, mol2

    def all_combinations_n_fittest(self, population, n):
        """
        Yields all pairings of the `n` fittest individuals.

        The pairings are yielded with in order of fitness, with most
        fit pair yielded first.

        Parameters
        ----------
        population : :class:`.GAPopulation`
            The population from which parents should be selected.

        n : :class:`int`
            The number of individuals used for making offspring.

        Yields
        ------
        :class:`tuple` of :class:`.MacroMolecule`
            The :class:`.MacroMolecule` instances which are to be
            crossed.

        """

        n_fittest = itertools.islice(self.fittest(population), n)
        for ind1, ind2 in itertools.combinations(n_fittest, 2):
            yield ind1, ind2

    def crossover_roulette(self, population, truncation=None):
        """
        Yields parents using roulette selection.

        In roulette selection the probability an individual is selected
        is given by its fitness. If the total fitness is the sum of all
        fitness values, the chance an individuals is selected is given
        by::

            p = individual fitness / total fitness,

        where ``p`` is the probability of selection [#]_.

        Parameters
        ----------
        population : :class:`.GAPopulation`
            The population from which parents are selected.

        truncation : :class:`int`, optional
            The number of least fit members which will never be
            yielded.

        Yields
        ------
        :class:`numpy.ndarray` of :class:`.MacroMolecule`
            The selected parent pair.

        References
        ----------
        .. [#] http://tinyurl.com/csc3djm

        """

        truncation = truncation if truncation is None else -truncation
        pop = sorted(population, reverse=True)[:truncation]

        while True:
            total = sum(ind.fitness for ind in pop)
            weights = [ind.fitness / total for ind in pop]
            yield np.random.choice(pop, 2, False, weights)

    def crossover_deterministic_sampling(self, population,
                                         truncation=None):
        """
        Yields parents according to determnistic sampling.

        In determnistic sampling the mean fitness value of the
        population is calculated, ``<f>``. For each individual a
        normalized fitness value is then calculated via::

            fn = f / <f>

        where ``fn`` is the normalized fitness value and ``f`` is the
        original fitness value.

        Deterministic sampling then creates a temporary, parent
        population of the same size as the original population. An
        individual is guaranteed to be placed into the parent
        population ``n`` times, where ``n`` is the integer part of
        ``fn``. Any remaining slots are given to individuals with the
        largest decimal values.

        Parents are then randomly selected from the parent population.

        Parameters
        ----------
        population : :class:`.GAPopulation`
            The population from which parents are selected.

        truncation : :class:`int`, optional
            The number of least fit members which will never be
            yielded.

        Yields
        ------
        :class:`numpy.ndarray` of :class:`.MacroMolecule`
            The selected parent pair.

        """

        truncation = truncation if truncation is None else -truncation
        pop = sorted(population, reverse=True)[:truncation]

        mean = np.mean([ind.fitness for ind in pop])

        parent_pop = []
        for ind in pop:
            fn = ind.fitness / mean
            if int(fn) < 1 and len(parent_pop) >= len(pop):
                break

            fn = 1 if int(fn) < 1 else int(fn)
            for x in range(fn):
                parent_pop.append(ind)

        while True:
            yield np.random.choice(parent_pop, 2, False)
