class FitnessCalculator:
    """
    Calculates fitness values of molecules.

    """

    def get_fitness(self, mol):
        """
        Return the fitness value of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule whose fitness value should be calculated.

        Returns
        -------
        :class:`object`
            The fitness value of `mol`.

        """

        return self._cache_result(self._get_fitness, mol)

    def _get_fitness(self, mol):
        """
        Return the fitness value of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule whose fitness should be calculated.

        Returns
        -------
        :class:`object`
            The fitness value of `mol`.

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method and needs to be implemented in a
            subclass.

        """

        raise NotImplementedError()


