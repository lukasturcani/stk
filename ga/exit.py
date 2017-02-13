"""
Defines functions which cut short the GA if a criterion is reached.

These functions are useful when debugging and testing the convergence
of the GA. For example they can tell the GA to stop when a certain
molecule has been found.

The functions are defined as methods in the Exit() class.

Extending MMEA: Adding Exit methods.
--------------------------------------
The only requirement is that the methods take `population` as their
first parameter (barring `self`, `cls` etc.) and return ``True`` if the 
exit criterion has been satisfied and ``False``` otherwise.

As usual, if you need to define multiple functions, make sure any 
helper functions are private, ie their names start with a leading 
underscore.

"""

class Exit:
    def __init__(self, func_data):
        self.func_data = func_data
        
    def __call__(self, pop):
        func = getattr(self, self.func_data.name)
        return func(pop, **self.func_data.params)

    def mol_present(self, population, mol):
        """
        ``True`` if `mol` is in `population`.
        
        Parameters
        ----------
        population : Population
            The GA population.
        
        mol : MacroMolecule
            A molecule which if present in `population` causes the GA to
            stop.
        
        Returns
        -------
        bool
            ``True`` if `mol` in `population`, ``False`` otherwise.
            
        """
        
        if mol in population:
            return True
        return False
        
    def no_exit(self, population):
        """
        Returns ``False``.
        
        Useful when you never want the GA to exit prematurely. This is 
        used by default when no exit function is defined in the input
        file.
        
        Returns
        -------
        False : bool
        
        """
        
        return False