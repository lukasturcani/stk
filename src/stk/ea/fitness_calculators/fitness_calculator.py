class FitnessCalculator:
    """
    Calculates fitness values of molecules.

    """

    def get_fitness_value(self, molecule):
        """
        Return the fitness value of `molecule`.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule whose fitness value should be calculated.

        Returns
        -------
        :class:`object`
            The fitness value of `molecule`.

        """

        return self._cache_result(self._get_fitness, mol)

    def _get_fitness(self, mol):
        """
        Return the fitness value of `molecule`.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule whose fitness value should be calculated.

        Returns
        -------
        :class:`object`
            The fitness value of `molecule`.

        """

        raise NotImplementedError()
