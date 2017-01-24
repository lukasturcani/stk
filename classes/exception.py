"""
A module for defining MMEA's exceptions.

"""

import os
import traceback

class MolError(Exception):
    """
    A class for raising errors when using ``Molecule`` instances.
    
    There are a lot of reason why MMEA might receive an error. The 
    error can originate when rkdit is trying to assembly a molucue, 
    sanitize it or manipulate it in some othe way. Equally, 
    optimizations may go wrong for one reason or another and raise an
    error. This is good, as problems should be identified, not ignored.
    
    However, errors raised by rdkit are not as useful as one would hope. 
    The error gets raised but the user is not told which ``Molecule`` 
    instance was being manipulated. This makes replication of the error 
    difficult.    
    
    In order to address this, all errors should be caught and placed
    into a ``MolError`` along with the Molecule instance on which the 
    error occured. On initialization of a MolError an entry is made in 
    the file ``failures.txt``. This file is located in the ``output`` 
    directory. Writing the entry to the file means that a Molecule can 
    be easily rebuilt as its building blocks, topology and any other 
    parameters are written to this file.
    
    Every try/except statement in  MMEA should be something like:
    
        try:
            raise SomeExceptionByRdkit(...)
        
        except Exception as ex:
            MolError(ex, mol, 'error in init')
            
            # Do some stuff here. Reraise if you want or pass.
    
    Notice that though an error is not raised, it is recorded in the
    ``failures.txt`` file. This is because a MolError instance was
    initialized during error handling.         
    
    Attributes
    ----------
    ex : Exception
        The exception raised by some other part of MMEA such as rdkit,
        or MacroModel.
        
    mol : Molecule
        The Molecule on which the error was raised.
        
    notes : str
        Any additional comments about the exception. For example, where
        it is being raised.
        
    """
    
    def __init__(self, ex, mol, notes):        
        self.ex = ex
        self.notes = notes
        self.write_to_file(mol)
        print('\n\nMolError written to ``failures.txt``.\n\n')
        
    def write_to_file(self, mol):
        """
        Writes the exception and Molecule to ``failures.txt``.
        
        This method is run during initialization. This means that even
        if an exception is ignored during runtime it will be still
        recorded in the ``failures.txt`` file.        
        
        """        
        
        # If the ``output`` folder exists (such as when running a GA
        # run) place the ``failures.txt`` file in it. If the file does
        # not exist (like when using MMEA as a library) place the
        # ``failures.txt`` in the same folder as ``MMEA``.
        cwd = os.getcwd().split('output')[0]
        if 'output' in os.getcwd():
            name = os.path.join(cwd, 'output', 'failures.txt')
        else:
            name = os.path.join(cwd, 'failures.txt')

        with open(name, 'a') as f:
            f.write("{} - {}\n\n".format(type(self.ex).__name__, 
                                            self.ex))
            
            traceback.print_exc(file=f)
            
            f.write('\nnote = {}\n'.format(self.notes))

            f.write('file = {}\n'.format(mol.file))
            
            if hasattr(mol, 'building_blocks'):       
                f.write('building blocks = {}\n'.format(
                                                   mol.building_blocks))
    
                f.write('topology = {}\n'.format(mol.topology)) 
               
                f.write('topology_args = {}\n'.format(
                                                     mol.topology_args))    
            f.write('\n'+'='*240)                 
            f.write('\n\n\n')                             
                             
class PopulationSizeError(Exception):
    def __init__(self, msg):
        self.msg = msg