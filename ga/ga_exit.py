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

    def __call__(self, pop, progress):
        func = getattr(self, self.func_data.name)
        return func(pop, progress, **self.func_data.params)

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

    def mol_name_present(self, population, progress, mol_name):
        """
        Returns ``True`` if molecule with `name` in `population`.

        Parameters
        ----------
        population : :class:`.Population`
            The GA population.

        progress : :class:`.Population`
            Population where each subpopulation is a previous generation.

        mol_name : :class:`str`
            The name of a molecule.

        Returns
        -------
        :class:`bool`
            ``True`` if a molecule with `name` of  `mol_name` is found
            in `population`.

        """

        return any(x.name == mol_name for x in population)

    def fitness_plateau(self, population, progress, num, top_members=1):
        """
        Returns ``True`` if fitness function of the top x candidates do not
        change for n generations.

        Parameters
        ----------
        population : :class:`Population`
            The GA population.

        progress : :class:`Population`
            Population where each subpopulation is a previous generation.

        num : :class:`int`
            Number of generations for which no the x top members did not
            change. This number is defined by the user.

        top_members : :class:`int`
            Number of members that are going to be considered for the plateau
            analysis. This number needs to be smaller than the total
            population. The user can define this number, which defaults to 1.

        Returns
        -------
        bool
            ``True`` if the fitness function shows no improvement for n
            generations.

        """
        # Check that the GA has run for more than num generations
        if len(progress.populations) > num:
            gens = {frozenset(sorted(progress.populations[-x],reverse=True)[:top_members]) for x in range(1, num)}
            unique_gens = len(gens)
            if unique_gens == 1:
                return True
            else:
                return False

        return False

    def no_exit(self, population, progress):
        """
        Returns ``False``.

        Useful when you never want the GA to exit prematurely. This is
        used by default when no exit function is defined in the input
        file.

        progress : :class:`Population`
            population where each subpopulation is a previous generation.

        Returns
        -------
        False : bool

        """

        return False
