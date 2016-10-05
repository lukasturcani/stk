from .ga import FunctionData
import MMEA.classes

def func_line_parser(line):
    """
    Creates a ``FunctionData`` instance based on data in line.
    
    This function must be applied only to lines which hold information
    about functions to be used by MMEA and their parameters.
    
    For details on what such a line should look like see the ``GAInput``
    class docstring.
    
    Parameters
    ----------
    line : str
        A line wihtin the MEA input file which defines a function and 
        its parameters.
        
    Returns
    -------
    FunctionData
        A ``FunctionData`` instance representing the MMEA function and 
        its parameters a defined within `line`.
    
    """
    
    # Split the line into components. Each component is text separated
    # by whitespace. Compenents = ``words``.
    words = line.split()
    # The first 2 words are the keyword defining the ``GAInput`` 
    # attribute name and the name of the function/method to be used,
    # respectively. Any remaining words are the parameters and their
    # values. At this point each parameter name and value is held 
    # together in a string in the form "param_name=param_val1" or 
    # param_name=param_val1,param_val2". All such parameter name and
    # value pairs are placed together in a list called ``params``.
    kw, name, *params = words
    # Create a dictionary which will hold the parameters name as the key
    # and the parameters value as the corresponding value.    
    param_dict = {}
    # Go through each parameter name-value pair and get each separately
    # by splitting at the ``=`` symbol. Then split any parameter values
    # joined by ``,`` using a split on this symbol. If only one 
    # parameter value was present the length of the list after the ``,``
    # split will be 1. In this case make just this one value the value 
    # in ``param_dict`` corresponding to the key. If more than one 
    # parameter value was present add all of them as the value in 
    # ``param_dict`` by making the value in ``param_dict`` a list 
    # holding all the values.    
    for param in params:
        p_name, p_vals = param.split("=")
        p_vals = p_vals.split(",")

        # If the parameter is a number, convert it to a float.
        for i, val in enumerate(p_vals):
            try:
                p_vals[i] = eval(val)
            except:
                pass

        if p_name == 'topologies':
            # Convert the strings of names of topologies into the actual
            # topology class objects. It is the classes themselves which
            # are used by MMEA not the class names as strings. So make
            # this conversion to facilitate that.            
            for i, topology in enumerate(p_vals):
                p_vals[i] = getattr(MMEA.classes.topology, p_vals[i])
            
        if len(p_vals) == 1 and isinstance(p_vals[0], str):
            # Sometimes file paths may include a space. Within the
            # input file the space in a path should be changed to 
            # ``~!~``. Here it is remade into a space.
            p_vals = p_vals[0].replace('~!~', ' ')                

        elif len(p_vals) == 1 and p_name != 'topologies':
            p_vals = p_vals[0]

        param_dict[p_name] = p_vals
        
    return FunctionData(name, **param_dict)
    
class GAInput:
    """
    A class for concisely holding information from MMEA's input file.
    
    A description of the input file follows, see also the User's guide.

    The input file consists of a sequence of lines. Each line defines a 
    variable or a function used by MMEA. If the line defines a function 
    used by MMEA the same line must also define any parameters necessary
    to use the function.
    
    If the line is empty or the first character is ``#`` it is skipped.
    This may be convenient if you wish to organize the input file into
    sections or add comments.
    
    Each non-empty line starts with a keyword. Each keyword corresponds
    to the name of one of the attributes defined in the ``Attributes``
    section of this docstring. For keywords which define a simple value
    such as ``pop_size`` or ``num_generations`` they are simply followed
    by the a space and the desired value.
    
    For lines where the keyword indicates a function is to being defined
    the syntax is as follows:
        
        keyword function_or_method_name param1_name=param1_value 
                param2_name=param2_value1,param2_value2
        
    NOTE: There would be no newline in the input file. The input file
          must be on a single line in its entirety. It was placed into
          2 in this docstring to account for MMEA style guidelines.
          
    Key points from the line example are:
        > Parameter names are followed by a ``=`` with NO WHITESPACE.
        > The ``=`` after the parameter name is followed by the value of
          the parameters with NO WHITESPACE.
        > If the parameter of a function expects an iterable such as a
          list or tuple, the values in that list or tuple are placed one
          after the other after the ``=``. Each value is separated from 
          the other by a ``,`` and NO WHITESPACE.
        > parameter names and value definitions are placed on the same
          line and separated by a whitespace.
          
    In some cases the parameter which needs to be provided will be a
    path. In cases where the path involves a space replace any space
    with ``~!~``. When the input file is read by this class it is remade
    into a space. For example if the path to MacroModel needs to be
    provided as a parameter:
    
        macromodel_path=C:\Program~!~Files\Schrodinger2016-2
    
    Notice the space in ``Program Files`` was replaced.
    
    Attributes
    ----------
    input_file : str
        The full path of the MMEA input file.
        
    pop_size : int
        The number of members in a MMEA generation.
        
    num_generations : int
        The number of generations formed by MMEA.
        
    num_mutations: int
        The number of successful mutations per generation.
        
    num_matings: int
        The number of successful matings per generation.
        
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
        
    mating_func : FunctionData
        The ``Mating`` class method used to mate ``MacroMolecule`` 
        instances to generate offspring. Must correspond to a method
        defined within the ``Mating`` class.
    
    mutation_func : FunctionData
        The ``Mutation`` class method used to mutate ``MacroMolecule`` 
        instances to generate mutants. Must correspond to a method 
        defined within the ``Mutation`` class.
        
    opt_func : FunctionData
        The function from the ``optimization.py`` module to be used for
        optimizing ``MacroMolecule`` instances.

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
        
        # If the input file did not specify the number of matings or
        # mutations it is assumed that none are wanted.
        if not hasattr(self, 'num_matings'):
            self.num_matings = 0
        
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
        # call the function which which exctracts data from function-
        # defining lines. If the keyword is not recognized, raise a 
        # ``NameError``.        
        with open(self.input_file, 'r') as input_file:
            for raw_line in input_file:
                line = raw_line.split()
                if not line or "#" == line[0]:
                    continue
                
                kw = line[0]
                
                if '_func' in kw:
                    func_data = func_line_parser(raw_line)
                    if 'mutation' in kw:
                        mutation_funcs = getattr(self, kw, [])
                        mutation_funcs.append(func_data)
                        func_data = mutation_funcs
                    setattr(self, kw, func_data)
                                        
                elif kw in {'pop_size', 'num_generations', 
                            'num_mutations', 'num_matings', 
                            'mutation_weights'}:
                    setattr(self, kw, eval(line[1]))
                
                else:
                    raise NameError(("Line does not define a valid"
                                                        " keyword."))
                    
                    