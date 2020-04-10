"""
Fitness Normalizer
==================

.. toctree::
    :maxdepth: 2

    Divide By Mean <stk.ea.fitness_normalizers.divide_by_mean>
    Multiply <stk.ea.fitness_normalizers.multiply>
    Null Fitness Normalizer <stk.ea.fitness_normalizers.null>
    Replace Fitness <stk.ea.fitness_normalizers.replace_fitness>
    Shift Up <stk.ea.fitness_normalizers.shift_up>
    Sum <stk.ea.fitness_normalizers.sum>

"""


class FitnessNormalizer:
    """
    Abstract base class for fitness normalizers.

    A fitness normalizer takes an :class:`tuple` of
    :class:`.MoleculeRecord` instances and yields new
    :class:`.MoleculeRecord` instances, with normalized fitness values.
    The primary benefit of a normalizer vs a
    :class:`.FitnessCalculator` is that a :class:`.FitnessNormalizer`
    has access to all members in the population when it is calculating
    the normalized fitness value, whereas a :class:`.FitnessCalculator`
    does not.

    """

    # Keep empty __init__() to override ugly default docstring.
    def __init__(self):
        """"""
        return

    def normalize(self, population):
        """
        Normalize the fitness values in `population`.

        Parameters
        ----------
        population : :class:`tuple` of :class:`.MoleculeRecord`
            The molecules which need to have their fitness values
            normalized.

        Yields
        ------
        :class:`.MoleculeRecord`
            A record with a normalized fitness value.

        """

        raise NotImplementedError()
