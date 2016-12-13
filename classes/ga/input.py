from types import ModuleType
import sys

from . import FunctionData
from ..topology import *
from ..population import Population
from .selection import Selection
from .crossover import Crossover
from .mutation import Mutation
from ...optimization import optimization
from ... import fitness

class GAInput:
    """
    A class for concisely holding information from MMEA's input file.
    
    A description of the input file follows, also see the User's guide.

    The input file consists of a sequence of commands. Each command 
    defines a variable or a function used by MMEA. If the command 
    defines a function used by MMEA it must also define any parameters 
    necessary to use the function. It does not have to define any 
    default initialized parameters, though it may if desired. A command 
    terminates with the ``$`` symbol. This means that 
    
        generational_select_func; 
        stochastic_sampling; 
        use_rank=True$
        
    and 
    
        generational_select_func; stochastic_sampling; use_rank=True$
        
    define the same command.
    
    If a line is empty or the first character is ``#`` it is skipped.
    This may be convenient if you wish to organize the input file into
    sections or add comments.
    
    Each non-empty line starts with a keyword. Each keyword corresponds
    to the name of one of the attributes defined in the ``Attributes``
    section of this docstring. For keywords which define a simple value
    such as ``num_generations`` they are simply followed by a ``=`` and 
    the desired value. For example,
        
        num_generations=25$
        
    would set the `num_generations` attribute of the ``GAInput`` 
    instance to 25. Notice there is no whitespace in this line. This is
    required.
    
    For commands where the keyword defines a function or method the 
    syntax is as follows:
        
        keyword; func_name; param1_name=param1_val; 
        param2_name=param2_val$
          
    Key points from the line example are:
        > Every unit is separated by a semicolon, ``;``, except the last
          which terminates with a ``$``.
        > Parameter names are followed by a ``=`` with NO WHITESPACE.
        > The ``=`` after the parameter name is followed by the value of
          the parameter with NO WHITESPACE.
    
    The ``func_name`` represents the name of a function or method which
    is being defined. For example:
    
        fitness_func; cage; target_cavity=5.7348; coeffs=[1,1,0,0,0]; 
        macromodel_path="/home/lukas/program_files/schrodinger2016-3"$

    NOTE: Not all parameters required by the ``cage`` function
          are defined.

    This command specifices that the ``cage()`` function defined within
    ``fitness.py`` is to be used as the fitness function. Notice that
    if the value passed to a parameter can be a list or a string.
    However, the type must be made explicit with either ``[]`` or quotes 
    for a string. Just like it would in a python script.
    
    Attributes
    ----------
    input_file : str
        The full path of the MMEA input file.
        
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
        instances for mutation from the current generation's population. 
        Must correspond to a method defined within the ``Selection`` 
        class.                        
        
    crossover_func : FunctionData
        The ``Crossover`` class method used to cross ``MacroMolecule`` 
        instances to generate offspring. Must correspond to a method
        defined within the ``Crossover`` class.
    
    mutation_func : FunctionData
        The ``Mutation`` class method used to mutate ``MacroMolecule`` 
        instances to generate mutants. Must correspond to a method 
        defined within the ``Mutation`` class.
        
    opt_func : FunctionData
        The function from the ``optimization.py`` module to be used for
        optimizing ``MacroMolecule`` instances.

    fitness_func : FunctionData
        The function from ``fitness.py`` to be used for calculating the
        fitness of ``MacroMolecule`` instances.
        
    """
    
    def __init__(self, input_file):
        """
        Initializes a ``GAInput`` instance.

        Parameters
        ----------
        input_file : str
            The full path of the MMEA input file.
            
        """
        
        self.input_file = input_file
        
        # Read the input file and extract its information.
        self._extract_data()
        
        # If the input file did not specify the number of crossovers or
        # mutations it is assumed that none are wanted.
        if not hasattr(self, 'num_crossovers'):
            self.num_crossovers = 0
        
        if not hasattr(self, 'num_mutations'):
            self.num_mutations = 0
            
        if not hasattr(self, 'mutation_weights'):
            self.mutation_weights = [1]
        
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
        
        Raises
        ------
        NameError
            If the keyword in the input file does not match any of the
            attribute names listed in the class level docstring.
        
        """
        
        # Open the input file and go through it line by line. If the 
        # keyword corresponds to a simple value just set it as the 
        # attribute and its value. If the keyword defines a function
        # call the function which extracts data from function defining 
        # lines. If the keyword is not recognized, raise a 
        # ``ValueError``.        
        with open(self.input_file, 'r') as input_file:
            
            # First remove all empty and comment lines.
            input_file = iter(line.strip() for line in input_file if 
                            not (line.isspace() or 
                                 line.strip()[0] == '#' or 
                                 line.strip() == ''))
                                 
            # Join up the file again and split across "$" to get full
            # commands.
            input_file = iter(line.strip() for line in 
                            " ".join(input_file).split("$") if
                            line != '')
                
            for raw_line in input_file:

                # Check if the keyword indicates a function defintion.
                kw, *_ = (word.strip() for word in raw_line.split(";"))
                if '_func' in kw:
                    func_data = self.line_data(raw_line)
                    if 'mutation' in kw:
                        mutation_funcs = getattr(self, kw, [])
                        mutation_funcs.append(func_data)
                        func_data = mutation_funcs
                    setattr(self, kw, func_data)
                    continue
            
                # Check if the keyword is a simple value. If it is, 
                # assign it to an attribute. If its not, raise a
                # ``ValueError``.
                try:
                    kw, val = raw_line.split("=")
                except Exception:
                    print(("\n\nERROR: Issue with the input file on the"
                    " following line (or its vicinity):\n\n"), raw_line,
                    "\n\n", sep="")
                    sys.exit()
        
                if kw in {'pop_size', 'num_generations', 'num_mutations', 
                          'num_crossovers', 'mutation_weights'}:
                    setattr(self, kw, eval(val))
                else:
                    raise ValueError(
                      "Line does not define a valid keyword.", raw_line)
                    
    @staticmethod
    def line_data(line):
        """
        Creates a ``FunctionData`` instance based on data in line.
        
        This function must be applied only to lines which hold 
        information about functions to be used by MMEA and their 
        parameters.
        
        For details on what such a line should look like see the 
        ``GAInput`` class docstring.
        
        Parameters
        ----------
        line : str
            A line wihtin the MEA input file which defines a function 
            and its parameters.
            
        Returns
        -------
        FunctionData
            A ``FunctionData`` instance representing the MMEA function 
            and its parameters as defined within `line`.
        
        """
        
        # Split the line into components. Each component is text 
        # separated by a semicolon. The components are essentially the
        # words on the line. The layout of a line is described in the 
        # class docstring above.
        kw, name, *params = (word.strip() for word in line.split(";"))
    
        # `param_dict` represents the parameters passed to the function
        # in `line` via the input file. It's a dictionary where the key 
        # is the name of a parameter defined in the input file  and the 
        # value is the corresponding value provided in the file. 
        param_dict = {}

        # Go through each parameter name-value pair in `line` and get 
        # each separately by splitting at the ``=`` symbol.   
        for param in params:
            try:
                p_name, p_vals = param.split("=")
            except ValueError:
                if param.count("=") > 1:
                    print(('\n\nERROR: Multiple "=" detected in the'
                    ' following line, did you forget a "$"?\n\n'), 
                    line, "\n\n", sep="")
                    sys.exit()
                
            except Exception:
                print(("\n\nERROR: Issue with the input file on the"
                " following line (or its vicinity):\n\n"), line, "\n\n",
                sep="")
                sys.exit()                        
            param_dict[p_name] = eval(p_vals)
            
        return FunctionData(name, **param_dict)
        
    def __repr__(self):
        return "\n\n".join("{} : {}".format(key, value) for key, value in
                         self.__dict__.items())
        
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
               'init_func' : iter(func for name, func in 
                                  Population.__dict__.items() if 
                                  name.startswith('init')),
               
               'generational_select_func' : iter(
                                 func for name, func in 
                                 Selection.__dict__.items() if 
                                 not name.startswith('crossover') and
                                 not name.startswith('_')),

               'parent_select_func' : iter(
                                 func for name, func in 
                                 Selection.__dict__.items() if 
                                 name.startswith('crossover')),
                                           
               'mutant_select_func' : iter(
                                 func for name, func in 
                                 Selection.__dict__.items() if 
                                 not name.startswith('crossover') and
                                 not name.startswith('_')),
                                           
               'crossover_func' : iter(func for name, func in 
                                    Crossover.__dict__.items() if 
                                    not name.startswith('_')),
               
               'mutation_func' : iter(func for name, func in 
                                    Mutation.__dict__.items() if 
                                    not name.startswith('_')),
                                      
               'opt_func' : iter(func for name, func in 
                                 optimization.__dict__.items() if
                                 not name.startswith('_') and 
                                 not isinstance(func, ModuleType) and
                                 'optimization' in func.__module__),
                                 
               'fitness_func' : iter(func for name, func in 
                                 fitness.__dict__.items() if
                                 not name.startswith('_') and 
                                 not isinstance(func, ModuleType) and 
                                 'fitness' in func.__module__)
               }
        
    
    def __init__(self, keyword):
        print('')
        for func in self.modules[keyword]:
            if hasattr(func, '__func__'):
                func = func.__func__
                
            print(func.__name__)
            print('-'*len(func.__name__))
            print(func.__doc__)
    
        
    
    
    
    
    
    
        