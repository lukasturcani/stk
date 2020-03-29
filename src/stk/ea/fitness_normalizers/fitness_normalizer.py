class FitnessNormalizer:
    """
    Abstract base class for fitness normalizers.

    A fitness normalizer takes an :class:`EAPopulation` and
    returns a new set of normalized fitness values for all
    molecules in it. The primary benefit of a normalizer vs a
    :class:`.FitnessCalculator` is that a :class:`.FitnessNormalizer`
    has access to all members in the population when it is calculating
    the new fitness value, whereas a :class:`.FitnessCalculator` does
    not.

    """

    def normalize(self, population):
        """
        Normalize the fitness values in `population`.

        Parameters
        ----------
        population : :class:`.EAPopulation`
            The molecules which need to have their fitness values
            normalized.

        Returns
        -------
        :class:`dict`
            Maps every molecule in `population` to its normalized
            fitness value.

        """

        # This method can be used to decorate _normalize in the future.
        return self._normalize(population)

    def _normalize(self, population):
        """
        Normalize fitness value in `population`.

        Parameters
        ----------
        population : :class:`.EAPopulation`
            The molecules which need to have their fitness values
            normalized.

        Returns
        -------
        :class:`dict`
            Maps every molecule in `population` to its normalized
            fitness value.

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method and needs to be implemented in a
            subclass.

        """

        raise NotImplementedError()
