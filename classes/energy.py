"""
Defines energy calculations via the ``Energy`` class.

Extending MMEA: Supporting more energy calculations.
----------------------------------------------------
It may well be the case that more ways to calculate energy of molcules
will need to be added. For example if support for using some external
software to do the calculations is needed.

In order to add a new function which calculates the energy, just add it 
as a method to the ``Energy`` class. There really aren't any 
requirements beyond this, just make sure the energy value is what the 
method returns.

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
    OUT: {('rdkit', {'forcefield' : 'uff'}) : 16.501}

IE this indicates that the `rdkit` method was called with the
`forcefield` argument set to 'uff'.    
    
Calling any method of the ``Energy`` class updates the dictionary 
automatically. When adding a new method to the class, no mechanism for
updating the dictionary needs to be provided. Writing the method within
the class is enough for it to update the dictionary when called.

Make use of this dictionary instead of running the same calculation 
repeatedly.

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

from .ga.cotainers import FunctionData

class EMethod:
    def __init__(self, func):
        self.func = func
        
    def __get__(self, obj, cls):
        
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
    
    def __new__(cls, cls_name, bases, cls_dict):
        
        # Find all the methods defined in the class `cls` and replace
        # them with an ``Emethod`` descriptor.
        for name, func in cls_dict.items():
            
            # Don't replace any private methods.
            if name.startswith('_') or not callable(func):
                continue
            
            # Replace method with descriptor.
            cls_dict[name] = EMethod(func)
        
        return type.__new__(cls, cls_name, bases, cls_dict)

class Energy(metaclass=EMeta):
    """
    Handles all things related to a ``Molecule``'s energy.

    An instance of this class will be placed in the `energy` attribute
    of a ``Molecule`` instance.    

    Attributes
    ----------
    molecule : Molecule    
        The energetic information held by an instance of ``Energy`` 
        concerns the molecule held in this attribute.
    
    values : dict
        The keys in the dict code for the function and parameters which
        were used to calculate the energy. The values are the energies.

    """

    def __init__(self, molecule):
        self.molecule = molecule
        self.values = {}
        
    def formation(self, key, products, 
                  building_blocks=None, force_e_calc=False,
                  nonkey_params=['force_e_calc']):
        """
        Calculates the formation energy.
        
        The formation energy is calculated under the assumption that the
        molecule in `self.molecule` is composed of the molecules in 
        `building_blocks` and that during formation molecules in 
        `products` are formed in addition to `self.molecule`.     
        
        Parameters
        ----------
        key : tuple
            The first member of the tuple is a string holding the name 
            of a method used to calculate energies. For exmaple 'rdkit'
            or 'macromodel'. The remaning elements in the tuple are the
            parameters that the user wishes to pass to the function.
            
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
        
        Modifies
        --------
        self.values : dict
            Adds an entry to this dictionary. The key is a tuple of the
            form ('formation', key[0], key[1]). The value is the 
            calculated formation energy.            
        
        Returns
        -------
        float
            The formation energy. Note that this value is also stored
            in the dictionary `self.values`.
        
        """

        func_name, *params = key
        
        # Recalculate energies if requested.
        if force_e_calc:
            for _, mol in products:
                getattr(mol.energy, func_name)(*params)
        
        e_products = 0
        for n, mol in products:
            if (func_name, params[0]) not in mol.energy.values.keys():
                getattr(mol.energy, func_name)(*params)
            
            e_products += n * mol.energy.values[(func_name, params[0])]
        
        eng = self.pseudoformation(key, building_blocks, force_e_calc) 
        eng -= e_products       
        return eng

    def pseudoformation(self, key, 
                        building_blocks=None,  force_e_calc=False, 
                        nonkey_params=['force_e_calc']):
        """
        Calculates the formation energy, sans other products.

        This is the formation energy if the energy of the other products
        of the reaction is not taken into account.
        
        Parameters
        ----------
        key : tuple
            The first member of the tuple is a string holding the name 
            of a method used to calculate energies. For exmaple 'rdkit'
            or 'macromodel'. The remaning elements in the tuple are the
            parameters that the user wishes to pass to the function.

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
        
        Modifies
        --------
        self.values : dict
            Adds an entry to this dictionary. The key is a tuple of the
            form ('formation', key[0], key[1]). The value is the 
            calculated formation energy sans products (pseudoformation).            
        
        Returns
        -------
        float
            The pseudoformation energy. Note that this value is also 
            stored in the dictionary `self.values`.
        
        """
        
        if building_blocks is None:
            building_blocks = ((n, mol) for mol, n in 
                              self.molecule.topology.bb_counter.items())
        
        func_name, *params = key
        
        # Recalculate energies if requested.
        if force_e_calc:
            for _, mol in building_blocks:
                getattr(mol.energy, func_name)(*params)
        
            getattr(self, func_name)(*params)

        # Calculate the energy of building blocks and products using the
        # chosen force field, if it has not been found already.
        e_reactants = 0
        for n, mol in building_blocks:
            if (func_name, params[0]) not in mol.energy.values.keys():
                getattr(mol.energy, func_name)(*params)
            
            e_reactants += n * mol.energy.values[(func_name, params[0])]
        
        e_products = (self.values[(func_name, params[0])] if 
                    (func_name, params[0]) in self.values.keys() else
                    getattr(self, func_name)(*params))

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
            ff = ac.MMFFGetMoleculeForceField(self.molecule.prist_mol,
                  ac.MMFFGetMoleculeProperties(self.molecule.prist_mol))

        eng = ff.CalcEnergy()        
        return eng
        
    def macromodel(self, forcefield, macromodel_path, 
                         nonkey_params=['macromodel_path']):
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
        key = call_sig(func, self, *args, **kwargs)
        
        # Update the `values` dictionary with the results of the 
        # calculation.
        obj.values.update({key:result})
        # Return the result.        
        return result

    # Make sure that the decorated function behaves like a typical 
    # method and return.
    return MethodType(inner, obj)

def call_sig(func, *args, **kwargs):
    
    fsig = sig(func)
    # Get a dictioanary of all the supplied parameters.
    bound = dict(fsig.bind_partial(*args, **kwargs).arguments)
    # Get a dictionary of all the default initialized parameters.
    default = {key : value.default for key,value in 
               dict(fsig.parameters).items() if key not in bound.keys()}
               
    # Combine the two sets of parameters and get rid of the `self` 
    # parameter.
    bound.update(default)
    bound.pop('self')
    # Remove any parameters that should not form key, listed in the 
    # `nonkey_params` parameter.
    if 'nonkey_params' in bound.keys():
        for key in bound['nonkey_params']:
            bound.pop(key)
        bound.pop('nonkey_params')
    
    # Return an FunctionData object representing the function and chosen
    # parameters.
    return FunctionData(func.__name__, **bound)
 


        