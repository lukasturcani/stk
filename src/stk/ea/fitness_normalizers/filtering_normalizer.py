
class _FilteringNormalizer(FitnessNormalizer):
    """
    Implements some of the :class:`.FitnessNormalizer` interface.

    """

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

        """

        def filter_fn(mol):
            return self._filter(population, mol)

        normalized = population.get_fitness_values()
        filtered = filter(filter_fn, population)
        # dict(normalized) is a copy of the initial dict, ensuring that
        # the dict used by _get_normalized_values does not change.
        normalized.update(
            self._get_normalized_values(filtered, dict(normalized))
        )
        return normalized

    def _get_normalized_values(self, filtered, fitness_values):
        """
        Yield normalized `filtered` fitness values.

        Parameters
        ----------
        filtered : :class:`iterable` of :class:`.Molecule`
            The molecules which passed the filter.

        fitness_values : :class:`dict`
            A mapping from every molecule in `filtered` to its
            fitness value.

        Yields
        ------
        :class:`tuple`
            Holds the :class:`.Molecule` as its first element and its
            normalized fitness value as its second.

        """

        raise NotImplementedError()


