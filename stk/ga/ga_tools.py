"""
Defines the :class:`.GATools`.

"""


class GATools:
    """
    Stores objects which carry out GA operations on populations.

    Attributes
    ----------
    selection: :class:`.Selection`
        The :class:`.Selection` instance which performs selections of a
        population's members.

    crossover : :class:`.Crossover`
        The :class:`.Crossover` instance which crosses a population's
        members.

    mutation : :class:`.Mutation`
        The :class:`.Mutation` instance which mutates a population's
        members

    normalization : :class:`.Normalization`
        A :class:`.Normalization` instance which rescales or normalizes
        the fitness values in the population.

    fitness : :class:`.FunctionData`, :class:`function`
        Holds the name and arguments of the function in :mod:`.fitness`
        to be used for calculating the fitness of
        :class:`.MacroMolecule` instances in the population.

        Alternatively, it can be a function which takes 1
        argument: the molecule whose fitness is to be calculated.

    exit : :class:`.Exit`
        An :class:`.Exit` object which checks if the population
        satisfied the exit criterion to stop the GA early.

    input : :class:`.GAInput`
        The :class:`.GAInput` object holding data gathered from the
        input file.

    """

    __slots__ = ['selection', 'crossover', 'mutation',
                 'normalization', 'fitness', 'exit', 'input']

    def __init__(self,
                 selection,
                 crossover,
                 mutation,
                 normalization,
                 fitness,
                 exit_,
                 ga_input):
        """
        Initializes a :class:`GATools` instance.

        Parameters
        ----------
        selection: :class:`.Selection`
            The :class:`.Selection` instance which performs selections
            of a population's members.

        crossover : :class:`.Crossover`
            The :class:`.Crossover` instance which crosses a
            population's members.

        mutation : :class:`.Mutation`
            The :class:`.Mutation` instance which mutates a
            population's members

        normalization : :class:`.Normalization`
            A :class:`.Normalization` instance which rescales or
            normalizes the fitness values in the population.

        fitness : :class:`.FunctionData`
            Holds the name and arguments of the function in
            :mod:`.fitness` to be used for calculating the fitness of
            :class:`.MacroMolecule` instances in the population.

        exit : :class:`.Exit`
            An :class:`.Exit` object which checks if the population
            satisfied the exit criterion to stop the GA early.

        input : :class:`.GAInput`
            The :class:`.GAInput` object holding data gathered from the
            input file.

        """

        self.selection = selection
        self.crossover = crossover
        self.mutation = mutation
        self.normalization = normalization
        self.fitness = fitness
        self.input = ga_input
        self.exit = exit_

    @classmethod
    def init_empty(cls):
        """
        Initializer which sets all attributes to ``None``.

        """

        return cls(None, None, None, None, None, None, None)
