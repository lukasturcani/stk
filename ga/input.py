"""
Defines classes which deal with input.

"""

from types import ModuleType
import sys
import re
from inspect import getmro
# Unused but may be used in input files. So needs to be present here as
# eval() is run parts of the input file.
import numpy as np

from . import fitness
from .crossover import Crossover
from .ga_tools import GATools
from .selection import Selection
from .mutation import Mutation
from .population import Population
from .normalization import Normalization
from .ga_exit import Exit

from ..convenience_tools import FunctionData
from ..molecular import topologies
from ..molecular.topologies import *
from ..molecular.molecules import *
from ..molecular import Energy
from ..molecular.optimization import optimization



class GAInput:
    """
    A class for concisely holding information from MMEA's input file.

    A description of the input file follows, also see the User's guide.

    The input file consists of a sequence of commands. Each command
    defines a variable or a function used by MMEA. If the command
    defines a function used by MMEA it must also define any parameters
    necessary to use the function. It does not have to define any
    default initialized parameters, though it may if desired. A command
    terminates at the start of the next command. Commands may be
    multiline, which means that

        generational_select_func;
        stochastic_sampling;
        use_rank=True

    and

        generational_select_func; stochastic_sampling; use_rank=True

    define the same command.

    If a line is empty or the first character is ``#`` it is skipped.
    This may be convenient if you wish to organize the input file into
    sections or add comments.

    Each non-empty line starts with a keyword. Each keyword corresponds
    to the name of one of the attributes defined in the ``Attributes``
    section of this docstring. For keywords which define a simple value
    such as ``num_generations`` they are simply followed by a ``=`` and
    the desired value. For example,

        num_generations=25

    would set the `num_generations` attribute of the ``GAInput``
    instance to 25. Notice there is no whitespace in this line. This is
    required.

    For commands where the keyword defines a function or method the
    syntax is as follows:

        keyword; func_name; param1_name=param1_val;
        param2_name=param2_val

    Key points from the line example are:
        > Every unit is separated by a semicolon, ``;``, except the
          last.
        > Parameter names are followed by a ``=`` with NO WHITESPACE.
        > The ``=`` after the parameter name is followed by the value
          of the parameter with NO WHITESPACE.

    The ``func_name`` represents the name of a function or method which
    is being defined. For example:

        fitness_func; cage; target_cavity=5.7348; coeffs=[1,1,0,0,0];
        macromodel_path="/home/lukas/program_files/schrodinger2016-3"

    This command specifices that the ``cage()`` function (defined
    within ``fitness.py``) is to be used as the fitness function.
    Notice that if the value passed to a parameter can be a list or a
    string. However, the type must be made explicit with either ``[]``
    or quotes for a string. Just like it would in a python script.

    If a new keyword is added to MMEA it should be added into the list
    `keywords`.

    The input file also allows definition of variables.

        $$hello=1
        $$world='hi'
        $$var3=FunctionData('somefunc', param1=$hello, param2=$world)

    The variables can hold anything. For example, $$var3 above can
    be re-written:

        $$PART1=FunctionData
        $$PART2=('somefunc', param1=$hello, param2=$world)
        $$var3=$PART1$PART2

    Variables do not have to be defined before they are used,

        $$var3=FunctionData('somefunc', param1=$hello, param2=$world)
        $$hello=1
        $$world='hi'

    is as valid as the first example.

    They essentially just perform the most basic text replacement.

    The variables are defined in the context of the input file, not
    within Python or the source code.

    To define a variable use ``$$`` followed by the name of the
    variable. To use a variable use `$` followed by the name of the
    variable.

    Class attributes
    ----------------
    _keywords : list
        Holds all valid keywords used by MMEA. Used to give users
        useful error messages.

    Attributes
    ----------
    _input_file : str
        The full path of the MMEA input file.

    pop_size : int
        The size of the population.

    num_generations : int
        The number of generations formed by MMEA.

    num_mutations: int
        The number of successful mutations per generation.

    num_crossovers: int
        The number of successful crossovers per generation.

    init_func : FunctionData
        The ``Population`` method used for initialization. This must
        correspond to a ``Population`` class initializer.

    generational_select_func : FunctionData
        The ``Selection`` class method used to select members of the
        next generation. Must correspond to a method defined within the
        ``Selection`` class.

    parent_select_func : FunctionData
        The ``Selection`` class method used to select parents from the
        current generation's population. Must correspond to a method
        defined within the ``Selection`` class.

    mutant_select_func : FunctionData
        The ``Selection`` class method used to select ``MacroMolecule``
        instances for mutation from the current generation's
        population. Must correspond to a method defined within the
        ``Selection`` class.

    crossover_func : FunctionData
        The ``Crossover`` class method used to cross ``MacroMolecule``
        instances to generate offspring. Must correspond to a method
        defined within the ``Crossover`` class.

    mutation_func : list of FunctionData instances
        The ``Mutation`` class methods used to mutate ``MacroMolecule``
        instances are held here. This is a list as multiple
        mutation functions can be used during the GAs run. The
        FunctionData instances mut correspond to a methods defined
        within the ``Mutation`` class.

    opt_func : FunctionData
        The function from the ``optimization.py`` module to be used for
        optimizing ``MacroMolecule`` instances.

    fitness_func : FunctionData
        The function from ``fitness.py`` to be used for calculating the
        fitness of ``MacroMolecule`` instances.

    mutation_weights : array-like
        The probability that each function in `mutation_func` will be
        selected each time a mutation operation is carried out. The
        order of the probabilities corresponds to the order of the
        mutation functions in `mutation_func`.

    normalization_func : list of FunctionData instances
        A list of functions which rescale or normalize the population's
        fitness values. The order reflects the order in which they are
        applied each generation.

    comparison_pops : list of strings
        A list of the full paths to pop_dump files which are to be
        compared. Only needed when using the `-c` option.

    databases : list of strings
        A list which holds the paths to any number JSON files. These
        files must hold the JSON represenatations of Population
        instances. All the molecules in the Populations are loaded
        into memory for the duration of the GA run. This means not all
        molecules have to be remade and optimized and have their
        fitness value recalculated.

    _content : str
        The string holding the data of the input file.

    _commands : list of str
        A list of all the commands in the input file.

    _variables : dict
        Maps all variables defined in the input file to their values.

    """

    _keywords = ['num_generations', 'num_mutations', 'num_crossovers',
                'init_func', 'generational_select_func', 'pop_size',
                'parent_select_func', 'mutant_select_func',
                'mutation_func', 'opt_func', 'mutation_weights',
                'crossover_func', 'fitness_func', 'normalization_func',
                'exit_func', 'comparison_pops', 'databases', r'\$\$']

    def __init__(self, input_file):
        """
        Initializes a ``GAInput`` instance.

        Parameters
        ----------
        input_file : str
            The full path of the MMEA input file.

        """


        self._variables= {}
        self._commands = []
        self._input_file = input_file
        with open(input_file, 'r') as inp:
            self._content = inp.read()

        # Replace any variables defined in the input file with the
        # corresponding values and any other pre-processing.
        self._process_input_file()
        # Read the input file and extract its information.
        self._extract_data()

        # If the input file did not specify some values, default
        # initialize them.
        if not hasattr(self, 'num_crossovers'):
            self.num_crossovers = 0

        if not hasattr(self, 'num_mutations'):
            self.num_mutations = 0

        if not hasattr(self, 'mutation_weights'):
            self.mutation_weights = [1]

        if not hasattr(self, 'normalization_func'):
            self.normalization_func = []

        if not hasattr(self, 'exit_func'):
            self.exit_func = FunctionData('no_exit')

        if not hasattr(self, 'databases'):
            self.databases = []

    @staticmethod
    def _command_data(command):
        """
        Create a FunctionData instance based on data in `command`.

        This function must be applied only to commands which hold
        information about functions to be used by MMEA and their
        parameters.

        For details on what such a command should look like see the
        ``GAInput`` class docstring.

        Parameters
        ----------
        command : str
            A command within the MEA input file which defines a
            function and its parameters.

        Returns
        -------
        FunctionData
            A ``FunctionData`` instance representing the MMEA function
            and its parameters as defined within `command`.

        """

        # Split the command into components. Each component is text
        # separated by a semicolon. The components are essentially the
        # words in the command. The layout of a command is described in
        # the class docstring above.
        kw, name, *params = (word.strip() for
                                        word in command.split(";"))

        # `param_dict` represents the parameters passed to the function
        # in `command` via the input file. It's a dictionary where the
        # key is the name of a parameter defined in the input file and
        # the value is the corresponding value provided in the file.
        param_dict = {}

        # Go through each parameter name-value pair in `command` and
        # get each separately by splitting at the ``=`` symbol.
        for param in params:
            p_name, p_vals = param.split("=", 1)
            param_dict[p_name] = eval(p_vals)

        return FunctionData(name, **param_dict)

    def crosser(self):
        """
        Returns a Crossover instance loaded with data from input file.

        Returns
        -------
        Crossover
            A Crossover instance which has all of the crossover
            related data held in the GAInput instance.

        """

        return Crossover(self.crossover_func, self.num_crossovers)

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

        return Exit(self.exit_func)

    def _extract_data(self):
        """
        Parses the input file and uses it to create attributes.

        Modifies
        --------
        self : GAInput
            Adds most of the attributes listed in the class docstring
            to the instance.

        Returns
        -------
        None : NoneType

        """

        # Go through each command. If the keyword corresponds to a
        # simple value just set it as the attribute and its value. If
        # the keyword defines a function call the function which
        # extracts data from function defining commands. If the keyword
        # is not recognized, raise a ``ValueError``.

        for command in self._commands:
            try:
                # Check if the keyword indicates a function
                # defintion.
                kw, *_ = (word.strip() for word in
                                            command.split(";"))
                if '_func' in kw:
                    func_data = self._command_data(command)

                    # In the case of mutation and normalization
                    # functions, place the extracted function into
                    # a list and then place that list as the
                    # attribute value.
                    if 'mutation' in kw or 'normalization' in kw:
                        funcs = getattr(self, kw, [])
                        funcs.append(func_data)
                        func_data = funcs

                    setattr(self, kw, func_data)
                    continue

                kw, val = command.split("=", 1)
                setattr(self, kw, eval(val))
            except:
                print(("\n\n\nERROR: Something is wrong with the"
                       " following command or in its vicinity.\n\n"),
                        command, sep="")
                sys.exit()

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
                       self.opt_func, self.fitness_func,
                       self.exiter(), self)

    def mutator(self):
        """
        Returns a Mutation instance loaded with data from input file.

        Returns
        -------
        Mutation
            A Mutation instance which has all of the mutation related
            data held in the GAInput instance.

        """

        return Mutation(self.mutation_func,
                        self.num_mutations, self.mutation_weights)

    def normalizer(self):
        """
        Returns Normalization instance holding data from input file.

        Returns
        -------
        Normalization
            A Normalization instance which has all of the normalization
            related data held in the GAInput instance.

        """

        return Normalization(self.normalization_func)

    def _process_input_file(self):
        """
        Does preprocessing on the input file.

        This includes things like substituting variables defined in
        the input file for their values and removing comments.

        Returns
        -------
        None : NoneType

        """

        # First remove all empty and comment lines.
        self._content = " ".join(line.strip() for line in
                            self._content.split('\n') if not
                            (line.strip() == '' or
                             line.isspace() or
                             line.strip()[0] == '#'))

        # Join up the file again and split across `keywords` to get
        # full commands and variables.
        p =  "(" + "|".join(self._keywords) + ")"
        p = re.compile(p)
        lines = [line for line in re.split(p, self._content) if line]
        keywords = lines[::2]
        content = lines[1::2]
        lines = [keyword + c for keyword, c in zip(keywords, content)]

        # Separate the commands and variables.
        for line in lines:
            if line.startswith('$$'):
                name, val = line[1:].split('=', 1)
                self._variables[name] = val
            else:
                self._commands.append(line)

        # Switch variables for their values.
        self._variable_sub()

    def selector(self):
        """
        Returns a Selection instance loaded with data from input file.

        Returns
        -------
        Selection
            A Selection instance which has all of the selection
            related data held in the GAInput instance.

        """

        return Selection(self.generational_select_func,
                         self.parent_select_func,
                         self.mutant_select_func)

    def _variable_sub(self, depth=0):
        """
        Substitutes variables for values.

        Parameters
        ----------
        depth : int (default = 0)
            This function is recursive. Depth keeps track of the
            recursion depth.

        Modifies
        --------
        _commands : list of strings
            The strings are changed so that all variables are replaced
            for the corresponding values.

        Raises
        ------
        RecursionError
            If `depth` exceeds 500. Possible when two variables form
            a cycle. Ie reference each other.

        Returns
        -------
        None : NoneType

        """

        content = "\n".join(self._commands)
        new_content = str(content)

        # Swith variables for values.
        for var in self._variables:
            new_content = new_content.replace(var,
                                              self._variables[var])

        # Remake the `_commands`.
        self._commands = new_content.split('\n')

        # If the new content is equal to the old content, all variables
        # have been substituted.
        if new_content == content:
            return

        # If not, and the recursion limit has been reached, exit.
        elif depth == 500:
            print(("\n\n\nRecursion limit reached while trying"
                   " to substitute variables. Two or more variables"
                   " are likely referencing each other.\n\n"))
            sys.exit()

        # If not, run the function again and try to substitute and
        # remaining nested variables.
        else:
            depth += 1
            return self._variable_sub(depth)

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
               'init_func' : (func for name, func in
                              Population.__dict__.items() if
                              name.startswith('init')),

               'generational_select_func' : (
                                 func for name, func in
                                 Selection.__dict__.items() if
                                 not name.startswith('crossover') and
                                 not name.startswith('_')),

               'parent_select_func' : (
                                 func for name, func in
                                 Selection.__dict__.items() if
                                 name.startswith('crossover')),

               'mutant_select_func' : (
                                 func for name, func in
                                 Selection.__dict__.items() if
                                 not name.startswith('crossover') and
                                 not name.startswith('_')),

               'crossover_func' : (func for name, func in
                                   Crossover.__dict__.items() if
                                   not name.startswith('_')),

               'mutation_func' : (func for name, func in
                                  Mutation.__dict__.items() if
                                  not name.startswith('_')),

               'opt_func' : (func for name, func in
                             optimization.__dict__.items() if
                             not name.startswith('_') and
                             not isinstance(func, ModuleType) and
                             'optimization' in func.__module__),

               'fitness_func' : (func for name, func in
                                 fitness.__dict__.items() if
                                 not name.startswith('_') and
                                 not isinstance(func, ModuleType) and
                                 'fitness' in func.__module__),

               'normalization_func' :  (func for name, func in
                                        Normalization.__dict__.items()
                                        if not name.startswith('_')),

                'energy' : (getattr(Energy, name) for name, func in
                                    Energy.__dict__.items() if not
                                    name.startswith('_')),

                'topologies' : (cls for name, cls in
                              topologies.__dict__.items() if
                              not name.startswith('_') and
                              not isinstance(cls, ModuleType) and
                              topologies.base.Topology in getmro(cls)),

                'exit_func' : (func for name, func in
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
