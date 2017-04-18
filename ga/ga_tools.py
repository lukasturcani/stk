"""
Defines the GATools class.

"""


class GATools:
    """
    Stores objects which carry out GA operations on populations.

    Instances of this class are held by ``Population`` instances in
    their `ga_tools` attribute. All the instances which carry out GA
    operations on the population are held conveniently together in an
    instance of this class. The optimization function used by the GA
    is also stored here as is the fitness function.

    Attributes
    ----------
    selection: Selection
        The ``Selection`` instance which performes selections of a
        population's members.

    crossover : Crossover
        The ``Crossover`` instance which crosses a population's members.

    mutation : Mutation
        The ``Mutation`` instance which mutates a population's members

    normalization : Normalization
        A ``Normalization`` instance which rescales or normalizes the
        fitness values of the population.

    optimization : FuncionData
        Holds the name of the function in ``optimization.py`` to be
        used for ``MacroMolecule`` optimizations and any additional
        parameters the function may require.

    fitness : FunctionData
        Holds the name of the function in ``fitness.py`` to be used for
        calculating the fitness of ``MacroMolecule`` instances. It also
        holds any additional paramters the function may require.

    exit : Exit
        An exit object which checks if the population achieved the
        exit criterion to stop the GA before all generations have been
        made.

    parallel : bool
        If ``True`` optimizations and fitness calculations are done on
        each molecule in parallel.

    input : GAInput
        The GAInput object holding data gathered from the input file.

    """

    __slots__ = ['selection', 'crossover', 'mutation', 'normalization',
                'optimization', 'fitness', 'exit', 'parallel', 'input']

    def __init__(self, selection, crossover,
                    mutation, normalization, optimization,
                         fitness, exit_, parallel, ga_input):
        self.selection = selection
        self.crossover = crossover
        self.mutation = mutation
        self.normalization = normalization
        self.optimization = optimization
        self.fitness = fitness
        self.input = ga_input
        self.parallel = parallel
        self.exit = exit_

    @classmethod
    def init_empty(cls):
        return cls(None, None, None, None,
                   None, None, None, None, None)
