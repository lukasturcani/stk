"""
Defines selection functions via the ``Selection`` class.

Extending MMEA: Adding selection functions
------------------------------------------
If a new selection operation is to be added to MMEA it should be added
as a method in the ``Selection`` class defined in this module. The only
requirements are that the first argument is ``population`` (excluding
any ``self`` or ``cls`` arguments) and that the method is decorated
wit the ``mutation()``, ``crossover()`` and ``generational()``
decorators to highlight which selection it is to be used for. The
decorators can be applied in sequence for functional which are suitable
for more than one kind of selection.

The naming requirement of ``population`` exists to help users identify
which arguments are handled automatically by MMEA and which they need
to defined in the input file. The convention is that if the mutation
function takes an argument called  ``macro_mol`` it does not have to be
specified in the input file.

All selection functios should be defined as generators, which yield an
member of ``population``. In the case of crossover selection functions,
they should yield a tuple of selected parents. Generational selection
function should yield each molecule at most once.

If the selection function does not fit neatly into a single function
make sure that any helper functions are private, ie that their names
start with a leading underscore.

"""

import itertools
import numpy as np


def mutation(func):
    """
    Identifies a function as suitable for mutation selections.

    Parameters
    ----------
    func : function
        A ``Selection`` class method which can be used for selecting
        molecules for mutation.

    Returns
    -------
    func : function
        The function received as an argument. The attribute
        ``mutation`` is added and set to ``True``.

    """

    func.mutation = True
    return func


def crossover(func):
    """
    Identifies a function as suitable for crossover selections.

    Parameters
    ----------
    func : function
        A ``Selection`` class method which can be used for selecting
        molecules for crossover.

    Returns
    -------
    func : function
        The function received as an argument. The attribute
        ``crossover`` is added and set to ``True``.

    """

    func.crossover = True
    return func


def generational(func):
    """
    Identifies a function as suitable for generational selections.

    Parameters
    ----------
    func : function
        A ``Selection`` class method which can be used for selecting
        molecules for the next generation.

    Returns
    -------
    func : function
        The function received as an argument. The attribute
        ``generational`` is added and set to ``True``.

    """

    func.generational = True
    return func


class Selection:
    """
    A class for handling all types of selection in the GA.

    Whenever a population needs to have some of its members selected
    for the creation of a parent pool or the next generation it
    delegates this task to an instance of this class. The population
    has this instance stored in its `ga_tools.selection` attribute.

    Each instance of this class supports being called. What a calling
    an instance does depends on the arguments the instance was
    initialized with and what arguments were supplied during the call.
    In all cases the call implements returns a generator. This
    generator yields members of a ``Population`` instance in accordance
    to some selection algorithm.

    Initialization of this class takes the names of methods defined in
    this class (as strings) and saves them into the instance's
    `generational`, `crossover` and `mutation` attributes. These
    attributes should therefore always hold the names of methods that
    are to be used for the given purpose - such as generational
    selection. The selection algorithms should be written as methods
    within this class.

    When an instance of this class is called it requires a
    ``Population`` instance and a string to be provided as arguments.
    The ``Population`` instance is the population which is to have some
    of its members selected. Consider the following code:

        >>> pop.select('generational')

    Here ``pop`` is a ``Population`` instance with a `ga_tools`
    attribute, which holds an initialized ``Selection`` instance.

    The method `select` invoked in the code automatically provides
    the instance ``pop`` to the ``Selection`` instance it contains. The
    function then calls the ``Selection`` instance. This means that
    each time `select` is called on a population, the ``Selection``
    instance will always act on the population it is held by. Different
    populations can use different selection algorithms by holding
    different ``Selection`` instances. Alternatively, they can perform
    the same selection algorithms by sharing a ``Selection`` instance.

    The string provided to the `select` method is passed all the way
    down to the ``Selection`` instance being called. This is the second
    argument a ``Selection`` instance requires when it is being called.
    This is the string `generational` in the code example above. The
    string should be the name of one of the attributes of the
    ``Selection`` class. This means that 'generational', 'crossover'
    and 'mutation' are valid strings at the time of this being written.
    If more types of selection are added to MMEA, an attribute named
    after that type of selection should be added to the ``Selection``
    class. If one wishes to invoke that type of selection on a
    population, `select` must be called with the name of that attribute
    as a parameter.

    There was a slight simplifcation in the paragraph describing
    initialization. When the ``Selection`` instance is initialzed it
    is not just with the names of the selection methods to be used. It
    provided with the names of the methods and any paramters that the
    methods will need to use. These parameters are packed into the
    ``FunctionData`` class. The ``FunctionData`` instance holding the
    method name and the appropriate parameter values is passed to
    the initializer of ``Selection``.

    Selection algorithms are to be implemented as generators. Selection
    algorithms which produce parents pools must yield a tuple of
    ``MacroMolecule`` instances. Selection algorithms should be grouped
    together by their expected use when written into the class body.

    Note that duplicates are not considered during selection. If a
    molecule is present in the population twice, its chance of
    selection does not double.

    Attributes
    ----------
    generational : FunctionData
        This holds the ``FunctionData`` object representing the
        selection function used for selecting the next generation of
        individuals along with any parameters the function may require.

    crossover : FunctionData
        Holds the ``FunctionData`` object representing the selection
        function used for selecting parents, along with any parameters
        the function may require.

    mutation : FunctionData
        Holds the ``FunctionData`` object representing the selection
        function used for selection individuals for mutation, along
        with any parameters the function may require.

    """

    def __init__(self, generational, crossover, mutation):
        self.generational = generational
        self.crossover = crossover
        self.mutation = mutation

    def __call__(self, population, type_):
        """
        Yields members from population.

        Parameters
        ----------
        population : Population
            The population from which members should be selected.

        type_ : str
            The name of a ``Selection`` attribute. The name corresponds
            to the type of selection that is desired: 'generational',
            'crossover' or 'mutation'.

        Yields
        -------
        MacroMolecule or tuple of MacroMolecules
            ``MacroMolecule`` instances are yielded unless `type_` is
            'crossover'. In this case a tuple of ``MacroMolecule``
            instances is yielded.

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

    @mutation
    @generational
    def fittest(self, population):
        """
        Yields members of the population, fittest first.

        Parameters
        ----------
        population : Population
            The population who's members need to be selected.

        Yields
        ------
        MacroMolecule
            The next fittest ``MacroMolecule`` instance held by the
            population.

        """

        for ind in sorted(population, reverse=True):
            yield ind

    @mutation
    @generational
    def roulette(self, population, elites=0,
                 truncation=None, duplicates=False):
        """
        Yields individuals using roulette selection.

        In roulette selection the probability an individual is selected
        is given by its fitness. If the total fitness is the sum of all
        fitness values, the chance an individuals is selected is given
        by

        p = individual fitness / total fitness,

        where p is the probability of selection [1].

        Parameters
        ----------
        population : Population
            The population from which individuals are to be selected.

        elites : int, optional
            The number of the fittest members which are guaranteed to
            be yielded first.

        truncation : int, optional
            The number of least fit members which will never be
            yielded.

        duplicates : bool, optional
            If ``True`` the same member can be yielded more than
            once.

        Yields
        ------
        MacroMolecule
            The next selected population member.

        References
        ----------
        [1] http://tinyurl.com/csc3djm

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

    @mutation
    def deterministic_sampling(self, population, truncation=None):
        """
        Yields individuals using deterministic sampling.

        This algorithm can only be used for selection of mutants. It
        will return duplicates. The results of ``fittest()`` are
        equivalent to algorithm without duplicates.

        In determnistic sampling the mean fitness value of the
        population is calculated, <f>. For each individual a
        normalized fitness value is then calculated via

            fn = f / <f>

        where fn is the normalized fitness value and f is the original
        fitness value. An individual will be selected n times, where
        n is the integer value of fn. After all individuals where
        n > 0 are yielded, all inidividuals are yielded again. Largest
        decimal part of fn first.

        Parameters
        ----------
        population : Population
            The population from which individuals are to be selected.

        truncation : int, optional
            The number of least fit members which will never be
            yielded.

        Yields
        ------
        MacroMolecule
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

    @mutation
    @generational
    def stochastic_sampling(self, population,
                            elites=0, truncation=None,
                            duplicates=False, use_rank=False):
        """
        Yields individuals via stochastic sampling.

        Each fitness value is used to calculate the normalized fitness
        of an individual

            fn = f / <f>

        where fn is the normalized fitness value, f is the original
        fitness value and <f> is the mean fitness of the population. If
        f is greater than 1 then the individual is guaranteed to be
        selected. If duplicates are allowed, individuals are guaranteed
        to be selected n times, where n is the integer part of fn.

        Once all the guarnteed individuals have been yielded, the
        remaining individuals are yielded via the roulette method. The
        weights in the roulette method are based on the decimal part of
        fn.

        Parameters
        ----------
        population : Population
            The population from which individuals are to be selected.

        elites : int, optional
            The number of the fittest members which are guaranteed to
            be yielded first.

        truncation : int, optional
            The number of least fit members which will never be
            yielded.

        duplicates : bool, optional
            If ``True`` the same member can be yielded more than
            once.

        use_rank : bool, optional
            When ``True`` the fitness value of an individual is
            calculated as, f = 1/rank.

        Yields
        ------
        MacroMolecule
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
        pop : list of MacroMolecules
            The molecules which are being selected.

        yielded : set of MacroMolecule instances
            Holds all previously yielded molecules.

        duplicates : bool
            Indicates whether a member can be yielded more than once.

        Yields
        ------
        MacroMolecule
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
        pop : list of MacroMolecules
            The molecules which are being selected.

        yielded : set of MacroMolecule instances
            Holds all previously yielded molecules.

        duplicates : bool
            Indicates whether a member can be yielded more than once.

        Yields
        ------
        MacroMolecule
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

    @crossover
    def all_combinations(self, population):
        """
        Yields every possible pairing of individuals from a population.

        This yields members regardless of which subpopulation they are
        in. Each pair is only returned once. This means if (1,2) is
        returned (2,1) will not be.

        Parameters
        ----------
        population : Population
            The population from which parents should be selected.

        Yields
        ------
        tuple of 2 MacroMolecule instances
            The ``MacroMolecule`` instances which together form a
            parent pair.

        """

        # Get all combinations of size 2.
        for mol1, mol2 in itertools.combinations(population, 2):
            yield mol1, mol2

    @crossover
    def all_combinations_n_fittest(self, population, n):
        """
        Yields all pairings of the `n` fittest individuals.

        The pairings are yielded with in order of fitness, with most
        fit pair yielded first.

        Parameters
        ----------
        population : Population
            The population from which parents should be selected.

        n : int
            The number of individuals used for making offspring.

        Yields
        ------
        tuple of 2 MacroMolecule instances
            The ``MacroMolecule`` instances which together form a parent
            pair.

        """

        n_fittest = itertools.islice(self.fittest(population), n)
        for ind1, ind2 in itertools.combinations(n_fittest, 2):
            yield ind1, ind2

    @crossover
    def crossover_roulette(self, population, truncation=None):
        """
        Yields parents using roulette selection.

        In roulette selection the probability an individual is selected
        is given by its fitness. If the total fitness is the sum of all
        fitness values, the chance an individuals is selected is given
        by

        p = individual fitness / total fitness,

        where p is the probability of selection [1].

        Parameters
        ----------
        population : Population
            The population from which parents are selected.

        truncation : int, optional
            The number of least fit members which will never be
            yielded.

        Yields
        ------
        numpy.ndarray of MacroMolecule instances
            The selected parent pair.

        References
        ----------
        [1] http://tinyurl.com/csc3djm

        """

        truncation = truncation if truncation is None else -truncation
        pop = sorted(population, reverse=True)[:truncation]

        while True:
            total = sum(ind.fitness for ind in pop)
            weights = [ind.fitness / total for ind in pop]
            yield np.random.choice(pop, 2, False, weights)

    @crossover
    def crossover_deterministic_sampling(self, population,
                                         truncation=None):
        """
        Yields parents according to determnistic sampling.

        In determnistic sampling the mean fitness value of the
        population is calculated, <f>. For each individual a
        normalized fitness value is then calculated via

            fn = f / <f>

        where fn is the normalized fitness value and f is the original
        fitness value.

        Deterministic sampling then creates a temporary, parent
        population of the same size as the original population. An
        individual is guaranteed to be placed into the parent
        population n times, where n is the integer part of fn. Any
        remaining slots are to individuals with the largest decimal
        values.

        Parents are then randomly selected from the parent population.

        Parameters
        ----------
        population : Population
            The population from which parents are selected.

        truncation : int, optional
            The number of least fit members which will never be
            yielded.

        Yields
        ------
        numpy.ndarray of MacroMolecule instances
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
