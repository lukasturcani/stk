"""
Defines energy calculations via the ``Energy`` class.

Extending MMEA: Supporting more energy calculations.
----------------------------------------------------
It may well be the case that more ways to calculate the energy of 
molcules will need to be added. For example, if support for using some 
external software to perform the calculations needs to be added.

In order to add a new function which calculates the energy, just add it 
as a method to the ``Energy`` class. There really aren't any 
requirements beyond this, just make sure the energy value is what the 
method returns. 

If the energy method does not fit neatly into a single method, feel free 
to split it up. Make sure sure that any helper methods are private, ie 
that their names start with a leading underscore. Only the main method
which the user will use should be public. 

To calculate the energy of some ``Molecule`` object the only thing that
is needed is calling one of the ``Energy`` methods:
    
    >>>  molecule.energy.rdkit('uff')
    OUT: 16.501
    
Here the `Energy.rdkit()` method was used as an example and `molecule` 
is just some  ``Molecule`` instance. (Each ``Molecule`` object has an 
``Energy`` instance in its `energy` attribute.)

Each ``Energy`` instance also has a `values` attribute which stores the 
result of previous calculations. Using the molecule from the previous 
example:

    >>>  molecule.energy.values
    OUT: {FunctionData('rdkit', forcefield=uff) : 16.501}

IE this indicates that the `rdkit` method was called with the
`forcefield` argument set to 'uff' and the result was 16.501.    
    
Calling any method of the ``Energy`` class updates the dictionary 
automatically. When adding a new method to the class, no mechanism for
updating the dictionary needs to be provided. Writing the method within
the class is enough for it to update the `values` dictionary when 
called.

Sometimes a function will require a parameter which does not affect the
outcome of the energy calculation. For example, to calculate the energy
using the ``macromodel`` program, the ``Energy.macromodel()`` function
can be used:
    
    molecule.energy.macromodel(forcefield=16, 
                               macromodel_path='path/to/macromodel/dir')
    
This function requires the number of a forcefield (16) and the directory
where ``macromodel`` is installed on the users computer 
('path/to/macromodel/dir'). However, the directory does not affect the
value of the calculated energy. When running:

    >>> molecule.energy.values
    
We want the output to be:
    
    OUT: {FunctionData('rdkit', forcefield=uff) : 16.501,
          FunctionData('macromolecule', forcefield=16): 200}

(Assuming we are still dealing with the same `molecule` instance from
the `rdkit` example, both calculated energies should be kept in the 
dictionary.)

However if we just define the ``Energy.macromodel()`` function within
the energy class, and then run:

    molecule.energy.macromodel(forcefield=16, 
                               macromodel_path='path/to/macromodel/dir')
                               
The output of 

    >>> molecule.energy.values

will be

    OUT: {FunctionData('rdkit', forcefield=uff) : 16.501,
          FunctionData('macromolecule', forcefield=16, 
                        macromodel_path='path/to/macromodel/dir'): 200}
                        
In order to make sure that the macromodel_path is excluded from the key,
decorate the function ``Energy.macromodel()`` with the ``exclude()`` 
decorator. For example:

    @exclude('macromodel_path')
    def macromodel(self, forcefield, macromodel_path):
        ...
        
Now the parameter `macromodel_path` will not form part of the key in
the `values` dictionary. If there were 2 parameters you wanted to 
exlude:

    @exclude('exclude1', 'exclude2')
    def energy_func(include1, include2, exclude1, exclude2):
        ...

and so on.

Exclusion of some parameters is doubly beneficial because it means less
typing needs to be done to make the key:

    >>>  key = FunctionData('macromodel', forcefield=16)
    >>>  molecule.energy.values[key]
    OUT: 200

rather than:

    >>>  key = FunctionData('macromodel', forcefield=16,
                             macromodel_path='path/to/macromodel/dir')
    >>>  molecule.energy.values[key]
    OUT: 200

Make sure to use the `values` dictionary instead of running the same 
calculation repeatedly.

The automatic updating of the dictionary is  achieved by the ``EMeta`` 
metaclass, ``EMethod`` descriptor and ``e_logger`` decorator. You do not 
need to worry about these, but information about how they work is 
provided in the docstring of ``EMeta``. 

"""

import os
import rdkit.Chem as chem
import rdkit.Chem.AllChem as ac
import subprocess as sp
import shutil
import uuid
from types import MethodType
from functools import wraps
from inspect import signature as sig

from .function_data import FunctionData

class EMethod:
    """
    A descriptor for methods of the ``Energy`` class.
    
    Attributes
    ----------
    func : function
        The method which the descriptor acts as a getter for.
    
    """
    
    def __init__(self, func):
        self.func = func
        
    def __get__(self, obj, cls):
        """
        Returns a modified `self.func`.
        
        Attributes
        ----------
        obj : object
            The object to which the method in `self.func` should be 
            bound.
            
        cls : object
            The class of `obj`.
        
        Returns
        -------
        BoundMethod
            A decorated version of the method held in `self.func`. The 
            difference is that when calling the method  now, it will 
            automatically update the `values` attribute of `obj`.
            
        self : EMethod
            If the method in `self.func` is called as a class attribute
            rather than an instance attribute, return the descriptor.
        
        """
        
        
        # If the Energy method is accessed as a class attribute return
        # the descriptor.
        if obj is None:
            return self
            
        # When trying to access the Energy method as an instance 
        # attribute returned a modified version of the method. The
        # method is modified so that after the method returns a value,
        # it is stored in the `values` dictionary of the Energy 
        # instance.
        return e_logger(self.func, obj)

class EMeta(type):
    """
    A metaclass for ``Energy``.    
    
    In conjuction with the EMethod descriptor and the e_logger decorator
    this function allows methods to automatically update the `values`
    attribute of their Energy instance without explicitly being told to
    do so.
    
    Basically this metaclass turns all methods of the Energy class into
    descriptors of the EMethod class. These descriptors return a 
    decorated version of the original method defined in the class. The 
    method is decorated with the ``e_logger()`` decorator. Calling this 
    decorated method makes it automatically update the `values` 
    dictionary of the object which used the method.
    
    """
    
    def __new__(cls, cls_name, bases, cls_dict):
        """
        Turns all the public methods in `cls` to EMethod descriptors.
        
        """
        
        # Find all the methods defined in the class `cls` and replace
        # them with an ``Emethod`` descriptor.
        for name, func in cls_dict.items():
            
            # Don't replace any private methods.
            if name.startswith('_') or not callable(func):
                continue
            
            # Replace method with descriptor.
            cls_dict[name] = EMethod(func)
        
        return type.__new__(cls, cls_name, bases, cls_dict)

def exclude(*args):
    """
    A decorator to add the `exclude` attribute to methods.

    Paremeters
    ----------
    args :  tuple of strings
        Holds the names parameters which are not to be used as part of
        the key identifying the energy calculation run.
        
    Returns
    -------
    function
        The function which has had the `exclude` attribute added. This
        is a list  holding the names of parameters of the function
        which are not used for identifying energy calculations.
    
    """

    
    def inner(func):
        func.exclude = args
        return func
        
    return inner

class Energy(metaclass=EMeta):
    """
    Handles all things related to a molecules energy.

    An instance of this class will be placed in the `energy` attribute
    of each ``Molecule`` instance.    

    Attributes
    ----------
    molecule : Molecule    
        The energetic information held by an instance of ``Energy`` 
        concerns the molecule held in this attribute.
    
    values : dict of FunctionData instances
        The keys in the dict code for the function and parameters which
        were used to calculate the energy. The values are the energies.

    """

    def __init__(self, molecule):
        self.molecule = molecule
        self.values = {}
    
    @exclude('force_e_calc')
    def formation(self, func, products, 
                  building_blocks=None, force_e_calc=False):
        """
        Calculates the formation energy.
        
        The formation energy is calculated under the assumption that the
        molecule in `self.molecule` is composed of the molecules in 
        `building_blocks` and that during formation molecules in 
        `products` are formed in addition to `self.molecule`.     
        
        Parameters
        ----------
        func : FunctionData
            A FunctionData object which describes the method of the 
            ``Energy`` class used to calculate the energies of the
            various molecules. For example:
                
                FunctionData('rdkit', forcefield='uff')
            
        products : tuple of form (float, Molecule)
            This tuple holds the molecules produced in addition to 
            `self.molecule`, when  a single `self.molecule` is made. The
            ``int`` represents the number made per `self.molecule`.

        building_blocks : tuple (default = None)
            This argument should be a tuple of the form 
            (float, Molecule). It holds the number of a given Molecule
            required to build a single molecule held in `self.molecule`.
            This argument can be omitted when the formation energy of a 
            MacroMolecule instance is being found, as they keep this 
            data stored elsewhere already.

        force_e_calc : bool (default = False)
            If the this is ``True`` then all building blocks, products
            and `self.molecule` will have their energies recalculated.
            Even if the energy values have already been found if the 
            chosen forcefield/method. If ``False`` the energy is only
            calculated if the value has not already been foud.
        
        Returns
        -------
        float
            The formation energy.
        
        """
        
        # Get the function used to calculate the energies.
        efunc = getattr(self, func.name)
        # Get the key of the function used to calculate the energies.
        fkey = func_key(efunc, **func.params)        
        
        # Recalculate energies if requested.
        if force_e_calc:
            for _, mol in products:
                getattr(mol.energy, func.name)(**func.params)
        
        # Get the total energy of the products.
        e_products = 0
        for n, mol in products:
            # If the energy has not been calculated already, calculate
            # it now.
            if fkey not in mol.energy.values.keys():
                getattr(mol.energy, func.name)(**func.params)
            
            e_products += n * mol.energy.values[fkey]
        
        eng = self.pseudoformation(func, building_blocks, force_e_calc) 
        eng -= e_products       
        return eng

    @exclude('force_e_calc')
    def pseudoformation(self, func, 
                        building_blocks=None,  force_e_calc=False):
        """
        Calculates the formation energy, sans other products.

        This is the formation energy if the energy of the other products
        of the reaction is not taken into account.
        
        Parameters
        ----------
        func : FunctionData
            A FunctionData object which describes the method of the 
            ``Energy`` class used to calculate the energies of the
            various molecules. For example:
                
                FunctionData('rdkit', forcefield='uff')

        building_blocks : tuple (default = None)
            This argument should be a tuple of the form 
            (float, Molecule). It holds the number of a given Molecule
            required to build a single molecule held in `self.molecule`.
            This argument can be omitted when the formation energy of a 
            MacroMolecule instance is being found, as they keep this 
            data stored elsewhere already.

        force_e_calc : bool (default = False)
            If the this is ``True`` then all building blocks, products
            and `self.molecule` will have their energies recalculated.
            Even if the energy values have already been found if the 
            chosen forcefield/method. If ``False`` the energy is only
            calculated if the value has not already been foud.        
        
        Returns
        -------
        float
            The pseudoformation energy.
        
        """

        # Get the function used to calculate the energies.
        efunc = getattr(self, func.name)
        # Get the key of the function used to calculate the energies.
        fkey = func_key(efunc, **func.params) 
        
        if building_blocks is None:
            building_blocks = ((n, mol) for mol, n in 
                              self.molecule.topology.bb_counter.items())
        
        # Recalculate energies if requested.
        if force_e_calc:
            for _, mol in building_blocks:
                getattr(mol.energy, func.name)(**func.params)
        
            getattr(self, func.name)(**func.params)

        # Sum the energy of building blocks under the chosen forcefield.
        e_reactants = 0
        for n, mol in building_blocks:
            # If the calculation has not been done already, perform it.
            if fkey not in mol.energy.values.keys():
                getattr(mol.energy, func.name)(**func.params)
            
            e_reactants += n * mol.energy.values[fkey]
        
        # Get the energy of `self.molecule`. The only product whose 
        # energy matters in pseudoformation.
        e_products = (self.values[fkey] if fkey in self.values.keys() 
                      else getattr(self, func.name)(**func.params))

        eng = e_reactants - e_products       
        return eng        

    def rdkit(self, forcefield):
        """
        Uses rdkit to calculate the energy of `self.molecule`.
        
        Parameters
        ----------
        forcefield : str
            The name of the forcefield to be used.
            
        Modifies
        --------
        self.values : dict
            Adds an entry to this dictionary. The key is a tuple of the
            form ('rdkit', `forcefield`). The value is the caculated 
            energy.
            
        Returns 
        -------
        float
            The calculated energy. Note that this value is also stored
            in the dictionary `self.values`.
        
        """
        
        if forcefield == 'uff':
            self.molecule.prist_mol.UpdatePropertyCache()
            ff = ac.UFFGetMoleculeForceField(self.molecule.prist_mol)
        if forcefield == 'mmff':
            chem.GetSSSR(self.molecule.prist_mol)
            self.molecule.prist_mol.UpdatePropertyCache()
            ff = ac.MMFFGetMoleculeForceField(self.molecule.prist_mol,
                  ac.MMFFGetMoleculeProperties(self.molecule.prist_mol))

        eng = ff.CalcEnergy()        
        return eng
     
    @exclude('macromodel_path')
    def macromodel(self, forcefield, macromodel_path):
        """
        Calculates the energy of `self.molecule` using macromodel.

        Note that this requires macromodel to be installed and have a 
        valid license.

        Parameters
        ----------
        forcefield : int
            The id number of the forcefield to be used by macromodel.  
            
        macromodel_path : str
            The full path of the ``Schrodinger`` suite within the user's 
            machine. For example, in a default Microsoft installation 
            the folder will probably be something like
            ``C:\Program Files\Schrodinger2016-2``.
 
        Modifies
        --------
        self.values : dict
            Adds an entry to this dictionary. The key is a tuple of the
            form ('macromodel', `forcefield`). The value is the 
            caculated energy.
            
        Returns 
        -------
        float
            The calculated energy. Note that this value is also stored
            in the dictionary `self.values`.
        
        """
        
        # To prevent conflicts when running this function in parallel,
        # a temporary copy of the molecular structure file is made and
        # used for macromodel calculations.
        
        # Unique file name is generated by inserting a random int into 
        # the prist_mol_file path.
        tmp_file = os.path.split(self.molecule.prist_mol_file)[-1]
        tmp_file, ext = os.path.splitext(tmp_file)
        r_int = uuid.uuid4().int
        tmp_file = tmp_file + str(r_int) + ext
        tmp_file = os.path.join(os.getcwd(), tmp_file)
        
        # Create the file.
        shutil.copy(self.molecule.prist_mol_file, tmp_file)
        
        file_root, ext = os.path.splitext(tmp_file)
        convrt_app = os.path.join(macromodel_path, 'utilities', 
                                                    'structconvert')
        convrt_cmd = [convrt_app, 
                     tmp_file, file_root+'.mae']
        sp.call(convrt_cmd, stdout=sp.PIPE, stderr=sp.PIPE)

        # Create an input file and run it.        
        input_script = (
        "{0}.mae\n"
        "{0}-out.maegz\n"
        " MMOD       0      1      0      0     0.0000     0.0000     "
        "0.0000     0.0000\n"
        " FFLD{1:8}      1      0      0     1.0000     0.0000     "
        "0.0000     0.0000\n"
        " BGIN       0      0      0      0     0.0000     0.0000     "
        "0.0000     0.0000\n"
        " READ      -1      0      0      0     0.0000     0.0000     "
        "0.0000     0.0000\n"
        " ELST      -1      0      0      0     0.0000     0.0000     "
        "0.0000     0.0000\n"
        " WRIT       0      0      0      0     0.0000     0.0000     "
        "0.0000     0.0000\n"
        " END       0      0      0      0     0.0000     0.0000     "
        "0.0000     0.0000\n\n"
        ).format(file_root, forcefield)
        
        with open(file_root+'.com', 'w') as f:
            f.write(input_script)
        
        cmd = [os.path.join(macromodel_path,'bmin'), 
               file_root, "-WAIT", "-LOCAL"]
        sp.call(cmd)
        
        # Read the .log file and return the energy.
        with open(file_root+'.log', 'r') as f:
            for line in f:
                if "                   Total Energy =" in line:
                    eng = float(line.split()[-2].replace("=", ""))
    
        
        # Clean up temporary files.
        for filename in os.listdir(os.path.split(tmp_file)[0]):
            if str(r_int) in filename:
                os.remove(filename)
        
        return eng
       
def e_logger(func, obj):
    
    @wraps(func)
    def inner(self, *args, **kwargs):
        
        # First get the result of the energy calculation.
        result = func(self, *args, **kwargs)
        
        # Next create FunctionData object to store the values of the
        # parameters used to to run that calculation.
        key = func_key(func, self, *args, **kwargs)
        
        # Update the `values` dictionary with the results of the 
        # calculation.
        obj.values.update({key:result})
        # Return the result.        
        return result

    # Make sure that the decorated function behaves like a typical 
    # method and return.
    return MethodType(inner, obj)

def func_key(func, *args, **kwargs):
    
    fsig = sig(func)
    # Get a dictionary of all the supplied parameters.
    bound = dict(fsig.bind_partial(*args, **kwargs).arguments)
    # Get a dictionary of all the default initialized parameters.
    default = {key : value.default for key,value in 
               dict(fsig.parameters).items() if key not in bound.keys()}
               
    # Combine the two sets of parameters and get rid of the `self` 
    # parameter, if present.
    bound.update(default)
    if 'self' in bound.keys():
        bound.pop('self')
    # Remove any parameters that should not form key, listed in the 
    # `exclude` attribute.
    if hasattr(func, 'exclude'):
        for key in func.exclude:
            bound.pop(key)
    
    # Return an FunctionData object representing the function and chosen
    # parameters.
    return FunctionData(func.__name__, **bound)
    
