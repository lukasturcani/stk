"""
Defines classes which deal with input.

"""

import psutil
import logging

from .ga_tools import GATools
from .selection import Selection
from .mutation import Mutation
from .crossover import Crossover
from .normalization import Normalization
from .ga_exit import Exit

from ..convenience_tools import FunctionData


class GAInput:
    """
    A class for reading a GA's input file.

    An ``stk`` input file is a Python script. The script must define a
    set of variables. Each variable defines a parameter or a function
    used by the GA. If the variable defines a function, it must also
    define any parameters necessary to use the function. It does
    not have to define any default initialized parameters, though it
    may if desired.

    Each variable name defined in the input file must corresond to an
    attribute of this class.

    Variables which define a function or method are supplied as
    dictionaries, for example:

    .. code-block:: python

        fitness_func = {'NAME' : 'func_name',
                        'param1_name' : param1_val,
                        'param2_name' : param2_val}

    Key points from the line example are:

        1. The variable name (e.g. ``fitness_func``) specifies the GA
           operation and matches the name of an attribute of
           :class:`GAInput`.
        2. The key ``'NAME'`` holds the name of the function which
           carries out the GA operation. Note the difference between
           the variable name specifying the generic role, i.e. a
           fitness function, versus the value for ``'NAME'`` holding
           the name of the function which carries out that role.
        3. Any parameters which the function needs are provided as key
           value pairs where the keys are strings holding the parameter
           name.

    Valid functions for a given role are described in the corresponding
    attribute of :class:`GAInput`.

    Some variables will need to define other parameters used by the GA
    such as constants or lists of constants,

    .. code-block:: python

        num_generations = 5
        databases = ['first/path', 'second/path']

    Attributes
    ----------
    processes : :class:`int`
        The number of parallel processes to create when running
        parallel operations. This does not have to be specified in the
        input file. If not specified, the number will default to the
        value returned by :func:`psutil.cpu_count`

    pop_size : :class:`int`
        The size of the population.

    num_generations : :class:`int`
        The number of generations the GA runs for.

    num_mutations: :class:`int`
        The number of mutation operations carried out per generation.

    num_crossovers: :class:`int`
        The number of crossover operations carried out per generation.

    init_func : :class:`dict`
        This must define the parameters of a :class:`.GAPopulation`
        initializer.

    generational_select_func : :class:`dict`
        This must define the parameters of a :class:`.Selection`
        method. This will be the function used to select members of the
        next generation.

    crossover_select_func : :class:`dict`
        This must define the parameters of a :class:`.Selection`
        method. This will be the function used to select molecules for
        crossover.

    mutation_select_func : :class:`dict`
        This must define the parameters of a :class:`.Selection`
        method. This will be the function used to select molecules for
        mutation.

    crossover_funcs : :class:`list` of :class:`dict`
        This :class:`list` holds a :class:`dict` for each crossover
        function which is to be used by the GA. The functions are
        defined as methods in :class:`.Crossover`

    mutation_funcs : :class:`list` of :class:`dict`
        This :class:`list` holds a :class:`dict` for each mutation
        function which is to be used by the GA. The functions are
        defined as methods in :class:`.Mutation`

    opt_func : :class:`dict`
        This must define the parameters of a function defined in
        :mod:`~stk.optimization`. This will be the function used to
        optimize the structures of the molecules.

    fitness_func : :class:`dict`
        This must define the parameters of a function defined in
        :mod:`.fitness`. This will be the function used to
        calculate the fitness of the molecules.

    mutation_weights : :class:`list` of :class:`int`
        The probability that each function in :attr:`mutation_funcs`
        will be selected each time a mutation operation is carried out.
        The order of the probabilities corresponds to the order of the
        mutation functions in :attr:`mutation_funcs`.

    normalization_funcs : :class:`list` of :class:`dict`
        This :class:`list` holds a :class:`dict` for each normalization
        function used to rescale or normalize the population's fitness
        values. The order reflects the order in which they are applied
        each generation. The functions are defined as methods in
        :class:`.Normalization`. This can be ommited from the input
        file which means that no normalization functions will be used.

    databases : :class:`list` of :class:`str`
        A :class:`list` which holds the paths to any number ``JSON``
        files. These files must hold the ``JSON`` representations of
        dumped :class:`.Population` instances. All the molecules in the
        populations are loaded into memory for the duration of the GA
        run. Molecules loaded in this way do not have to be built,
        re-optimized or have fitness values recalculated when
        encountered by the GA.

    progress_dump : :class:`bool`
        If ``True`` a ``.json`` :class:`.Population` dump is made at
        the end of the GA run called ``progress.json``. The population
        holds each generation made by the GA as a subpopulation. This
        can be ommitted from the input file, in which case it
        defaults to ``True``.

    database_dump : :class:`bool`
        If ``True`` a ``.json`` :class:`.Population` dump is made at
        the end of the GA run called ``database.json``. The population
        holds every molecule made by the GA during its run. This can
        be omitted from the input file, in which case it defaults to
        ``True``.

    plot_epp : :class:`bool` or :class:`str`
        If a string is provided, it should hold the base name for the
        EPP plots made by the GA. If ``False`` then no EPP plots are
        made. This can be omitted from the input file, in which case
        it defaults to ``'epp.png'``.

    logging_level : :class:`int`
        The logging level for logging messages to the screen. This can
        be omitted from the input file in which case it defulats to
        ``logging.DEBUG``.

    counters : :class:`bool`
        If ``True`` plots of which molecules get chosen by selection
        functions are made. This can be omitted from the input file,
        in which case it defaults to ``True``.

    pop_dumps : :class:`bool`
        If ``True`` each generation makes a ``.json`` dump file. This
        can be omitted from the input file, in which case it defaults
        to ``True``. It means that a progress and database dump is made
        at each generation as well.

    tar_output : :class:`bool`
        If ``True`` a copy of the output folder is compressed, tarred
        and placed in the output folder. This can be omitted from the
        input file, in which case it defaults to ``True``.

    progress_load : :class:`str`
        The path to a ``progress.json`` file generated by a previous
        GA run. This means the previous GA progress is restored and a
        restart GA run can be performed. Can be omitted from the input
        file in which case it dafaults to ``''``. This means that a
        new GA run is assumed.

    """

    def __init__(self, input_file):
        """
        Initializes a :class:`GAInput` instance.

        Parameters
        ----------
        input_file : :class:`str`
            The path to the GA input file.

        """

        with open(input_file, 'r') as inp:
            exec(inp.read(), vars(self))

        if (hasattr(self, 'plot_epp') and
           getattr(self, 'plot_epp') is True):
            self.plot_epp = 'epp.png'

        # If the input file did not specify some values, default
        # initialize them.
        if not hasattr(self, 'processes'):
            self.processes = psutil.cpu_count()

        if not hasattr(self, 'num_crossovers'):
            self.num_crossovers = 0

        if not hasattr(self, 'num_mutations'):
            self.num_mutations = 0

        if not hasattr(self, 'mutation_weights'):
            self.mutation_weights = None

        if not hasattr(self, 'crossover_weights'):
            self.crossover_weights = None

        if not hasattr(self, 'normalization_funcs'):
            self.normalization_funcs = []

        if not hasattr(self, 'exit_func'):
            self.exit_func = {'NAME': 'no_exit'}

        if not hasattr(self, 'databases'):
            self.databases = []

        if not hasattr(self, 'progress_dump'):
            self.progress_dump = True

        if not hasattr(self, 'database_dump'):
            self.database_dump = True

        if not hasattr(self, 'plot_epp'):
            self.plot_epp = 'epp.png'

        if not hasattr(self, 'logging_level'):
            self.logging_level = logging.DEBUG

        if not hasattr(self, 'counters'):
            self.counters = True

        if not hasattr(self, 'pop_dumps'):
            self.pop_dumps = True

        if not hasattr(self, 'tar_output'):
            self.tar_output = True

        if not hasattr(self, 'progress_load'):
            self.progress_load = ''

    def crosser(self):
        """
        Use input file data to build :class:`.Crossover` instance.

        Returns
        -------
        :class:`.Crossover`
            A :class:`.Crossover` instance which has all of the
            crossover related data held by the :class:`.GAInput`
            instance.

        """

        funcs = [FunctionData(x['NAME'],
                              **{k: v for k, v in x.items() if
                                 k != 'NAME'})
                 for x in self.crossover_funcs]

        return Crossover(funcs,
                         self.num_crossovers,
                         self.crossover_weights)

    def exiter(self):
        """
        Use input file data to build :class:`.Exit` instance.

        Returns
        -------
        :class:`.Exit`
            An :class:`.Exit` instance loaded with the exit function
            defined in the input file.

        """

        func_data = FunctionData(self.exit_func['NAME'],
                                 **{key: val for key, val in
                                    self.exit_func.items() if
                                    key != 'NAME'})
        return Exit(func_data)

    def fitnessor(self):
        """
        Use input file data to  extract fitness function information.

        Returns
        -------
        :class:`.FunctionData`
            A :class:`.FunctionData` object which represents the
            fitness function.

        """

        return FunctionData(self.fitness_func['NAME'],
                            **{key: val for key, val in
                               self.fitness_func.items() if
                               key != 'NAME'})

    def ga_tools(self):
        """
        Use input file data to build a :class:`.GATools` instance.

        Returns
        -------
        :class:`.GATools`
            A :class:`.GATools` instance built from the information
            defined in the input file.

        """

        return GATools(self.selector(),
                       self.crosser(),
                       self.mutator(),
                       self.normalizer(),
                       self.fitnessor(),
                       self.exiter(),
                       self)

    def initer(self):
        """
        Use input file to  extract initialization function information.

        Returns
        -------
        :class:`.FunctionData`
            A :class:`.FunctionData` object which represents the
            initialization function.

        """

        return FunctionData(self.init_func['NAME'],
                            **{key: val for key, val in
                               self.init_func.items() if
                               key != 'NAME'})

    def mutator(self):
        """
        Use input file data to build :class:`.Mutation` instance.

        Returns
        -------
        :class:`.Mutation`
            A :class:`.Mutation` instance which has all of the
            mutation related data held in the input file.

        """

        funcs = [FunctionData(x['NAME'],
                 **{k: v for k, v in x.items() if k != 'NAME'})
                 for x in self.mutation_funcs]

        return Mutation(funcs,
                        self.num_mutations,
                        self.mutation_weights)

    def normalizer(self):
        """
        Use input file data to build :class:`.Normalization` instance.

        Returns
        -------
        :class:`.Normalization`
            A :class:`.Normalization` instance which has all of the
            normalization related data held in the input file.

        """

        funcs = [FunctionData(x['NAME'],
                 **{k: v for k, v in x.items() if k != 'NAME'})
                 for x in self.normalization_funcs]
        return Normalization(funcs)

    def opter(self):
        """
        Use the input file to extract optimization function data.

        Returns
        -------
        :class:`.FunctionData`
            A :class:`.FunctionData` object which represents the
            optimization function defined in the input file.

        """

        return FunctionData(self.opt_func['NAME'],
                            **{key: val for key, val in
                               self.opt_func.items() if key != 'NAME'})

    def selector(self):
        """
        Use input file data to build :class:`.Selection` instance.

        Returns
        -------
        :class:`.Selection`
            A :class:`.Selection` instance which has all of the
            selection related data held in the input file.

        """

        gen = FunctionData(self.generational_select_func['NAME'],
                           **{key: val for key, val in
                              self.generational_select_func.items() if
                              key != 'NAME'})

        crossover = FunctionData(self.crossover_select_func['NAME'],
                                 **{key: val for key, val in
                                    self.crossover_select_func.items()
                                    if key != 'NAME'})

        mutation = FunctionData(self.mutation_select_func['NAME'],
                                **{key: val for key, val in
                                   self.mutation_select_func.items() if
                                   key != 'NAME'})

        return Selection(gen, crossover, mutation)

    def __repr__(self):
        return "\n\n".join("{} : {}".format(key, value) for key, value
                           in self.__dict__.items())

    def __str__(self):
        return repr(self)
