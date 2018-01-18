"""
Defines functions which cut short the GA if a condition is satisfied.

These functions are useful when debugging and testing the convergence
of the GA. For example they can tell the GA to stop when a certain
molecule has been found.

The functions are defined as methods in the :class:`Exit`.

Extending mtk: Adding exit functions.
-------------------------------------

The only requirement is that the methods take `population` as their
first argument (barring `self`, `cls` etc.) and `progress` as their
second. They must return ``True`` if the exit criterion has been
satisfied and ``False`` otherwise. The naming requirement of the first
two arguments exists to help users identify which arguments they need
to define in the input file and which are handeled automatically by
``mtk``. The arguments `population` and `progress` are automatically
provided when the GA is run. They hold a :class:`.Population` instance
holding the current generation and a :class:`.Population` instance
holding every previous generation in a subpopulation, respectively.

As usual, if you need to define multiple functions, make sure any
helper functions are private, i.e. their names start with a leading
underscore.

"""


class Exit:
    """
    Checks if the exit criterion for the GA has been satisfied.

    Attributes
    ----------
    func_data : :class:`.FunctionData`
        A :class:`.FunctionData` object which holds the name of the
        exit function to use and any additional parameters it may
        require.

    """

    def __init__(self, func_data):
        """
        Initializes an :class:`Exit` instance.

        Parameters
        ----------
        func_data : :class:`.FunctionData`
            A :class:`.FunctionData` object which holds the name of the
            exit function to use and any additional parameters it may
            require.

        """

        self.func_data = func_data

    def __call__(self, pop, progress):
        """
        Calls the chosen exit function.

        Parameters
        ----------
        pop : :class:`.Population`
            The population holding the current generation of molecules.

        progress : :class:`.Population`
            A population where every previous generation is a
            subpopulation.

        Returns
        -------
        :class:`bool`
            The value returned by giving `pop` and `progress` to the
            exit function defined in :attr:`func_data`

        """

        func = getattr(self, self.func_data.name)
        return func(pop, progress, **self.func_data.params)

    def mol_present(self, population, progress, mol):
        """
        ``True`` if `mol` is in `population`.

        Parameters
        ----------
        pop : :class:`.Population`
            The population holding the current generation of molecules.

        progress : :class:`.Population`
            A population where every previous generation is a
            subpopulation.

        mol : :class:`.MacroMolecule`
            A molecule which if present in `population` causes the GA
            to stop.

        Returns
        -------
        :class:`bool`
            ``True`` if `mol` in `population`, ``False`` otherwise.

        """

        if mol in population:
            return True
        return False

    def mol_name_present(self, population, progress, mol_name):
        """
        Looks for a molecule called `mol_name`.

        Parameters
        ----------
        pop : :class:`.Population`
            The population holding the current generation of molecules.

        progress : :class:`.Population`
            A population where every previous generation is a
            subpopulation.

        mol_name : :class:`str`
            The name of a molecule.

        Returns
        -------
        :class:`bool`
            ``True`` if a molecule with :attr:`.MacroMolecule.name` of
            `mol_name` is found in `population`. Else, ``False``.

        """

        return any(x.name == mol_name for x in population)

    def fitness_plateau(self,
                        population,
                        progress,
                        num_gens,
                        top_members=1):
        """
        Checks if the most fit molecules remain the same.

        Parameters
        ----------
        pop : :class:`.Population`
            The population holding the current generation of molecules.

        progress : :class:`.Population`
            A population where every previous generation is a
            subpopulation.

        num_gens : :class:`int`
            Number of generations in which the fittest molecules did
            not change.

        top_members : :class:`int`, optional
            The number of fittest molecules which are checked. This
            number needs to be smaller than the population size.

        Returns
        -------
        :class:`bool`
            ``True`` if the fitness function shows no improvement for n
            generations.

        """

        # Check that the GA has run for more than num_gens generations.
        if len(progress.populations) > num_gens:
            gens = {frozenset(sorted(progress.populations[-x],
                                     reverse=True)[:top_members])
                    for x in range(1, num_gens)}
            unique_gens = len(gens)
            if unique_gens == 1:
                return True
            else:
                return False

        return False

    def no_exit(self, population, progress):
        """
        Returns ``False``.

        Useful when you never want the GA to exit prematurely.

        Parameters
        ----------
        pop : :class:`.Population`
            The population holding the current generation of molecules.

        progress : :class:`.Population`
            A population where every previous generation is a
            subpopulation.

        Returns
        -------
        False : :class:`bool`

        """

        return False
