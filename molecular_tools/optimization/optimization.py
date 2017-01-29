"""
Defines optimization functions.

Extending MMEA: Adding optimization functions
---------------------------------------------
New optimization functions are added by writing them into this module.
The only requirements are that the first argument is ``macro_mol`` and
that the function returns the optimized instance of ``macro_mol``. The
frist requirement allows users to identify which arguments are handled
automatically by MMEA and which need to be defined in the input file.
The convention is that if the optimization function takes an argument 
called ``macro_mol`` the users does not have to specify that argument in 
the input file. The second requirment is necessary for reasons to do
with parallization. If it is not fulfilled, the molecule will not be
updated with the optimized version.

Optionally optimizations can be complicated. If the use of helper
functions is required make sure that they are private, ie that their
names begin with a leading underscore. In the event that the
optimization is so complex that it requires its own module, place the 
module in the to same folder as this file. Then import the optimization
function into this file. See the ``macromodel.py`` as an example. Make
sure that only the optimization functions are imported back into this
file, not any of the helper functions or classes. 

"""

import rdkit.Chem.AllChem as ac
import rdkit.Chem as chem
import multiprocessing as mp
from functools import partial

from .macromodel import (macromodel_opt, 
                         macromodel_cage_opt, macromodel_md_opt)

def _optimize_all(func_data, population):
    """
    Apply optimization function to all population members in parallel.

    Individual optimization functions defined within this module should 
    change the `optimized` attribute of ``MacroModel`` instances to 
    ``True``. They should also include the line
    
        if macro_mol.optimized:
            return None
    
    at the start. This prevents optimizing an already optimized
    structure again.
    
    If this function should be used, rather than its serial counterpart
    ``optimize_all_serial()``, the ``optimize_population`` method in the
    ``Population`` class must be told to use it.

    The parallel optimization creates cloned instances of the 
    population's members. It is these that are optimized. This means 
    that the structure files are changed but any instance attributes 
    are not.

    To deal with this, optimization functions should return the 
    ``MacroMolecule`` instance they optimize. The pool will return an 
    iterator of the values returned by the optimization functions. If 
    the returned values are the optimized macromolecules, the population 
    in the main thread can be updated with them.

    Parameters
    ----------
    func_data : FunctionData
        The ``FunctionData`` object which represents the chosen
        optimization function. This function should be defined within
        this module. The ``FunctionData`` object also holds any
        additional parameters the optimization function may need.
        
    population : Population
        The ``Population`` instance who's members must be optimized.
        
    Modifies
    --------
    MacroMolecule's structure files
        This function optimizes the structures of all the 
        ``MacroMolecule`` instances held in `population`. This means
        that their structure files are modified to their optimized 
        structures. However, because all instances are cloned
        any values in the original instance's attributes are unchanged.
        The file's contents are changed because they are written on the
        hard disk, which can be clones of the python interpreter.
    
    Returns
    -------
    iterator of MacroMolecule objects
        This iterator yields the ``MacroMolecule`` objects that have had
        their attributes changed as a result of the optimization. They
        are modified clones of the original population's macromolecules.
    
    """
    
    # Using the name of the function stored in `func_data` get the
    # function object from one of the functions defined within the 
    # module.
    func = globals()[func_data.name]
    # Provide the function with any additional paramters it may require.
    p_func = partial(func, **func_data.params)
    
    # Apply the function to every member of the population, in parallel.
    with mp.get_context('spawn').Pool() as pool:
        optimized = pool.map(p_func, population)
        # Make sure the cache is updated with the optimized versions.
        for member in optimized:
            member.update_cache()
        return optimized
    
def _optimize_all_serial(func_data, population):
    """
    Apply optimization function to all population members, serially.

    Individual optimization functions defined within this module should 
    change the `optimized` attribute of ``MacroModel`` instances to 
    ``True``. They should also include the line
    
        if macro_mol.optimized:
            return None
    
    at the start. This prevents optimizing an already optimized
    structure again.
    
    If this function should be used, rather than its parallel 
    counterpart ``optimize_all``, the ``optimize_population`` method in 
    the ``Population`` class must be told to use it.

    Parameters
    ----------
    func_data : FunctionData
        The ``FunctionData`` object which represents the chosen
        optimization function. This function should be defined within
        this module. The ``FunctionData`` object also holds any
        additional parameters the optimization function may need.
        
    population : Population
        The ``Population`` instance who's members must be optimized.
        
    Modifies
    --------
    MacroMolecule
        This function optimizes the structures of all the 
        ``MacroMolecule`` instances held in `population`. This means
        that their structure files are modified to their optimized 
        structures. However, only the content of these files is changed. 
        The value of the `file` attributes remain the same. 
    
    Returns
    -------
    iterator of MacroMolecule objects
        This is meant to mirror the output of the parallel counterpart.
        This allows the two functions to be interfaced in the same way.
    
    """

    # Using the name of the function stored in `func_data` get the
    # function object from one of the functions defined within the 
    # module.    
    func = globals()[func_data.name]
    # Provide the function with any additional paramters it may require.
    p_func = partial(func, **func_data.params)
    
    # Apply the function to every member of the population.    
    return iter(p_func(member) for member in population)    
    
def rdkit_optimization(macro_mol):
    """
    Optimizes the structure of the molecule using rdkit.
    
    Parameters
    ----------
    macro_mol : MacroMolecule
        The macromolecule who's structure should be optimized.
        
    Modifies
    --------
    macro_mol.mol
        The rdkit molecule held in this attribute has it's structure
        changed as a result of the optimization. This means the
        ``Conformer`` instance held by the rdkit molecule is changed.
    
    macro_mol.file's content
        The content of the structure file located at 
        `macro_mol.file`, is changed so that it holds the structure of 
        the optimized rdkit molecule.
    
    macro_mol.optimized
        After a successful optimization, this attribute is set to 
        ``True``.
    
    Returns
    -------
    macro_mol : MacroMolecule   
        The macromolecule that was passed as an argument and modified
        by the optimization function. Returned to accomodate 
        parallelization. See ``optimize_all`` function documentation for
        more details.
    
    """
    
    # If `macro_mol` is already optmized, return.
    if macro_mol.optimized:
        print('Skipping {0}.'.format(macro_mol.file))   
        return macro_mol
        
    # Sanitize then optimize the rdkit molecule.
    chem.SanitizeMol(macro_mol.mol)
    ac.MMFFOptimizeMolecule(macro_mol.mol)
    
    # Update the content of the structure file.
    macro_mol.write()
    
    macro_mol.optimized = True   
    return macro_mol
    
def do_not_optimize(macro_mol):
    """
    Skips the optimization step.
    
    This is very useful when debugging so you do not waste your time
    waiting for molecules to get optimized. Use this in the input file
    in place of an optimization function when necessary.
    
    Parameters
    ----------
    macro_mol : MacroMolecule
        A macromolecule which will not be optimized.
        
    Modifies
    --------
    macro_mol.optimized
        Set to ``True``.
    
    Returns
    -------
    MacroMolecule
        The macromolecule not getting optimized.
    
    """
    
    if macro_mol.optimized:
        print('Skipping', macro_mol.file)
        return macro_mol
    
    print('Optimizing', macro_mol.file)
    macro_mol.optimized = True   
    return macro_mol   
