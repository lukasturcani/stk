"""
Defines classes which deal with input.

"""

from types import ModuleType
from inspect import isclass
import logging

from . import fitness
from .crossover import Crossover
from .ga_tools import GATools
from .selection import Selection
from .mutation import Mutation
from .population import Population
from .normalization import Normalization
from .ga_exit import Exit

from .. import ga
from .. import molecular
from .. import convenience_tools

from ..convenience_tools import FunctionData
from ..molecular import Energy, topologies
from ..molecular.optimization import optimization


class GAInput:
    """
    A class for concisely holding information from MMEA's input file.

    An MMEA input file is a Python script. The script must define a set
    of variables. Each variable defines a parameter or a function used
    by MMEA. If the variable defines a function used by MMEA it must
    also define any parameters necessary to use the function. It does
    not have to define any default initialized parameters, though it
    may if desired.

    Each variable name must corresond to an attribute of this class.

    Variables which define a function or method are supplied as
    dictionaries, for example:

        fitness_func = {'NAME' : 'func_name',
                        'param1_name' : param1_val,
                        'param2_name' : param2_val}

    Key points from the line example are:
        > The variable name specifies the GA operation and is equal to
          the name of an attribute of this class.
        > The key 'NAME' holds the name of the function which carries
          out the GA operation.
        > Any parameters which the function needs are provided as key
          value pairs where the keys are strings holding the parameter
          name.

    Valid function names for a given variable can be found by using
    the -h option.

        python -m mmea -h fitness_func

    See also the User's guide.

    Some variables will need to define other parameters used by MMEA,
    such as constants,

        num_generations = 5

    or lists of constants

        databases = ['first/path', 'second/path']

    See the ``Attributes`` section to see what data each variable must
    hold.

    Attributes
    ----------
    pop_size : int
        The size of the population.

    num_generations : int
        The number of generations formed by MMEA.

    num_mutations: int
        The number of successful mutations per generation.

    num_crossovers: int
        The number of successful crossovers per generation.

    init_func : dict
        The key 'NAME' must hold the name of a ``Population`` method
        initializer method.

    generational_select_func : dict
        The key 'NAME' must hold the name of the ``Selection`` class
        method used to select members of the next generation.

    crossover_select_func : dict
        The key 'NAME' must hold the name of the ``Selection`` class
        method used to select parents from the current generation's
        population for crossover.

    mutation_select_func : dict
        The key 'NAME' must hold the name of the ``Selection`` class
        method used to select ``MacroMolecule`` instances for mutation
        from the current generation's population.

    crossover_funcs : list of dicts
        This list holds a dict for each crossover function which is to
        be used by the GA. The key 'NAME' in each dict must hold the
        name of a method of the ``Crossover`` class.

    mutation_funcs : list of dicts
        This list holds a dict for each mutation function which is to
        be used by the GA. The key 'NAME' in each dict must hold the
        name of a method of the ``Mutation`` class.

    opt_func : dict
        The key 'NAME' of the dict must hold the name of a fucntion
        defined in the ``optimization.py`` module. It is used for
        optimizing the structure of generated molecules.

    fitness_func : dict
        The key 'NAME' must hold the name of a function defined in
        ``fitness.py``. The function is used to calculate the fitness
        of generated molecules.

    mutation_weights : list of ints
        The probability that each function in `mutation_funcs` will be
        selected each time a mutation operation is carried out. The
        order of the probabilities corresponds to the order of the
        mutation functions in `mutation_funcs`.

    normalization_funcs : list of dicts
        A list of functions which rescale or normalize the population's
        fitness values. The order reflects the order in which they are
        applied each generation.

    databases : list of strings
        A list which holds the paths to any number JSON files. These
        files must hold the JSON representations of Population
        instances. All the molecules in the populations are loaded
        into memory for the duration of the GA run. This means not all
        molecules have to be remade and optimized or have their
        fitness value recalculated.

    parallel : bool (default = True)
        If ``True`` MMEA does optimization and fitness calculation of
        each molecule in parallel.

    progress_dump : bool (default = True)
        If ``True`` a .json Population dump is made at the end of the
        GA run called ``progress`.json`. The population holds each
        generation as a subpopulation.

    database_dump : bool (default = True)
        If ``True`` a .json Population dump is made at the end of the
        GA run called ``database.json``. The population holds every
        molecule made by the GA as a member.

    plot_epp : False or str (default = 'epp.png')
        If a string, then it should hold the name of the EPP plots.
        If ``False`` then no EPP plots are made.

    logging_level : int (default = logging.DEBUG)
        The logging level for logging messages to the screen.

    counters : bool (default = True)
        If ``True`` plots the which molecules get chosen by selection
        functions get made.

    pop_dumps : bool (default = True)
        If ``True`` each generation makes a .json dump file.

    tar_output : bool (default = True)
        If ``True`` a copy of the output folder is compressed, tarred
        and placed in the output folder.

    """

    def __init__(self, input_file):
        """
        Initializes a ``GAInput`` instance.

        Parameters
        ----------
        input_file : str
            The full path of the MMEA input file.

        """

        with open(input_file, 'r') as inp:
            mmea_namespace = {}
            mmea_namespace.update(vars(convenience_tools))
            mmea_namespace.update(vars(molecular))
            mmea_namespace.update(vars(ga))
            exec(inp.read(), mmea_namespace, vars(self))

        if (hasattr(self, 'plot_epp') and
           getattr(self, 'plot_epp') is True):
            self.plot_epp = 'epp.png'

        # If the input file did not specify some values, default
        # initialize them.
        if not hasattr(self, 'num_crossovers'):
            self.num_crossovers = 0

        if not hasattr(self, 'num_mutations'):
            self.num_mutations = 0

        if not hasattr(self, 'mutation_weights'):
            self.mutation_weights = [1]

        if not hasattr(self, 'crossover_weights'):
            self.crossover_weights = [1]

        if not hasattr(self, 'normalization_funcs'):
            self.normalization_funcs = []

        if not hasattr(self, 'exit_func'):
            self.exit_func = {'NAME': 'no_exit'}

        if not hasattr(self, 'databases'):
            self.databases = []

        if not hasattr(self, 'parallel'):
            self.parallel = True

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

    def crosser(self):
        """
        Returns a Crossover instance loaded with data from input file.

        Returns
        -------
        Crossover
            A Crossover instance which has all of the crossover
            related data held in the GAInput instance.

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
        Returns a Exit instance loaded with data from the input file.

        Returns
        -------
        Exit
            An Exit instance loaded with the exit function defined in
            the input file. If none was defined an exit function which
            always returns ``False`` is used.

        """

        func_data = FunctionData(self.exit_func['NAME'],
                                 **{key: val for key, val in
                                 self.exit_func.items() if
                                 key != 'NAME'})
        return Exit(func_data)

    def fitnessor(self):
        """
        Returns a FunctionData of fitness func in input file.

        Returns
        -------
        FunctionData
            A FunctionData object which represents the fitness
            function.

        """

        return FunctionData(self.fitness_func['NAME'],
                            **{key: val for key, val in
                            self.fitness_func.items() if
                            key != 'NAME'})

    def ga_tools(self):
        """
        Return a GATools instance loaded with data from the input file.

        Returns
        -------
        GATools
            A GATools instance which has all of the input data held in
            the GAInput instance.

        """

        return GATools(self.selector(), self.crosser(),
                       self.mutator(), self.normalizer(),
                       self.opter(), self.fitnessor(),
                       self.exiter(), self.parallel, self)

    def initer(self):
        """
        Returns a FunctionData object of init function in input file.

        Returns
        -------
        FunctionData
            A FunctionData object which represents the initialization
            function.

        """

        return FunctionData(self.init_func['NAME'],
                            **{key: val for key, val in
                            self.init_func.items() if
                            key != 'NAME'})

    def mutator(self):
        """
        Returns a Mutation instance loaded with data from input file.

        Returns
        -------
        Mutation
            A Mutation instance which has all of the mutation related
            data held in the GAInput instance.

        """

        funcs = [FunctionData(x['NAME'],
                 **{k: v for k, v in x.items() if k != 'NAME'})
                 for x in self.mutation_funcs]

        return Mutation(funcs,
                        self.num_mutations,
                        self.mutation_weights)

    def normalizer(self):
        """
        Returns Normalization instance holding data from input file.

        Returns
        -------
        Normalization
            A Normalization instance which has all of the normalization
            related data held in the GAInput instance.

        """

        funcs = [FunctionData(x['NAME'],
                 **{k: v for k, v in x.items() if k != 'NAME'})
                 for x in self.normalization_funcs]
        return Normalization(funcs)

    def opter(self):
        """
        Returns a FunctionData of optimization func in input file.

        Returns
        -------
        FunctionData
            A FunctionData object which represents the optimization
            function.

        """

        return FunctionData(self.opt_func['NAME'],
                            **{key: val for key, val in
                            self.opt_func.items() if
                            key != 'NAME'})

    def selector(self):
        """
        Returns a Selection instance loaded with data from input file.

        Returns
        -------
        Selection
            A Selection instance which has all of the selection
            related data held in the GAInput instance.

        """

        gen = FunctionData(self.generational_select_func['NAME'],
                           **{key: val for key, val in
                           self.generational_select_func.items() if
                           key != 'NAME'})

        crossover = FunctionData(self.crossover_select_func['NAME'],
                                 **{key: val for key, val in
                                 self.crossover_select_func.items() if
                                 key != 'NAME'})

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


class InputHelp:
    """
    A class which creates output when ``-h`` option is used as input.

    The ``-h`` option is used in the following way:

        python -m MMEA -h keyword

    Here ``keyword`` corresponds to one of the attributes of the
    ``GAInput`` class. The output when this command is used will be
    a list of all functions which can be used with that keyword and
    the corresponding documentation.

    Class attributes
    ----------------
    modules : dict
        Maps the name of the keyword to the object which holds the
        functions or methods that are to be used with that keyword.

    """

    modules = {
               'init_func': (func for name, func in
                             Population.__dict__.items() if
                             name.startswith('init')),

               'generational_select_func': (
                                 func for name, func in
                                 Selection.__dict__.items() if
                                 hasattr(func, 'generational')),

               'crossover_select_func': (
                                 func for name, func in
                                 Selection.__dict__.items() if
                                 hasattr(func, 'crossover')),

               'mutation_select_func': (
                                 func for name, func in
                                 Selection.__dict__.items() if
                                 hasattr(func, 'mutation')),

               'crossover_funcs': (func for name, func in
                                   Crossover.__dict__.items() if
                                   not name.startswith('_')),

               'mutation_funcs': (func for name, func in
                                  Mutation.__dict__.items() if
                                  not name.startswith('_')),

               'opt_func': (func for name, func in
                            optimization.__dict__.items() if
                            not name.startswith('_') and
                            not isinstance(func, ModuleType) and
                            'optimization' in func.__module__),

               'fitness_func': (func for name, func in
                                fitness.__dict__.items() if
                                not name.startswith('_') and
                                not isinstance(func, ModuleType) and
                                'fitness' in func.__module__),

               'normalization_funcs':  (func for name, func in
                                        Normalization.__dict__.items()
                                        if not name.startswith('_')),

               'energy': (getattr(Energy, name) for name, func in
                          Energy.__dict__.items() if not
                          name.startswith('_')),

               'topologies': (cls for name, cls in
                              topologies.__dict__.items() if
                              not name.startswith('_') and
                              isclass(cls) and
                              issubclass(cls, topologies.base.Topology)),

               'exit_func': (func for name, func in
                             Exit.__dict__.items() if not
                             name.startswith('_'))
               }

    def __init__(self, keyword):
        print('')
        for func in self.modules[keyword]:
            if hasattr(func, '__func__'):
                func = func.__func__

            print(func.__name__)
            print('-'*len(func.__name__))
            print(func.__doc__)
