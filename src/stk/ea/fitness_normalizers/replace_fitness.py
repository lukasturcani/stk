"""
Replace Fitness
===============

"""

from .fitness_normalizer import FitnessNormalizer


class ReplaceFitness(FitnessNormalizer):
    """
    Replaces fitness values of a certain value with a new value.

    Examples
    --------
    Replace all fitness values which are ``None`` with the half the
    minimum fitness value in the population

    .. code-block:: python

        import stk

        replacer = stk.ReplaceFitness(
            replacement_fn=lambda population:
                min(
                    f for _, f in population.get_fitness_values()
                    if f is not None
                ) / 2,
            filter=lambda population, mol:
                population.get_fitness_values()[mol] is None,
        )

    """

    def __init__(
        self,
        get_replacement,
        filter=lambda population, mol: True,
    ):
        """
        Initialize a :class:`.ReplaceFitness` instance.

        Parameters
        ----------
        get_replacement : :class:`callable`
            Takes a single parameter, the :class:`.Population` which
            needs to be normalized, before it is filtered, and
            returns an :class:`object` which is used as the new
            fitness value for all molecules which pass the
            `filter`.

        filter : :class:`callable`, optional
            Takes a two parameters, first is a :class:`.EAPopulation`
            and the second is a :class:`.Molecule`, and
            returns ``True`` or ``False``. Only molecules which
            return ``True`` will have fitness values replaced. By
            default, all molecules will have fitness values replaced.
            The :class:`.EAPopulation` on which :meth:`normalize` is
            called is passed as the first argument while the second
            argument will be passed every :class:`.Molecule` in it.

        """

        self._get_replacement = get_replacement
        self._filter = filter

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

        def filter_fn(mol):
            return self._filter(population, mol)

        replacement_value = self._replacement_fn(population)
        normalized = population.get_fitness_values()
        for mol in filter(filter_fn, population):
            normalized[mol] = replacement_value
        return normalized
