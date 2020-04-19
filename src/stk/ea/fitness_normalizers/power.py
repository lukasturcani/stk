"""
Power
=====

"""

import numpy as np

from .fitness_normalizer import FitnessNormalizer


class Power(FitnessNormalizer):
    """
    Raises fitness values to some power.

    Examples
    --------
    *Raising Fitness Values to a Power*

    Sometimes you might calculate a property for a molecule, where
    that property indicates a low fitness value. You can use
    :class:`.Power` to raise it to the power of -1 to get your
    final fitness value

    .. code-block:: python

        normalizer = stk.Power(-1)
        # Assuming that population holds molecule record instances
        # with the following fitness values: 1, 2, 3
        # normalized will hold fitness values of
        # 1, o.5, 1/3
        normalized = tuple(normalizer.normalize(population))


    *Raising Fitness Values by a Set of Powers*

    In this example, assume that each fitness value consists of a
    :class:`tuple` of numbers, each representing a different property
    of the molecule, and each contributing to the final fitness value.
    The properties can be anything, such as  energy, number of atoms
    or diameter.

    If your final fitness value depends on the combination of these
    properties, you will probably first want to scale them with
    :class:`.DivideByMean`. Once this is done, you may want to
    raise each property by some power. For example
    if you raise the value of one property by ``1`` and another
    by ``-1``, it means that when the value of property 1 is big,
    the fitness value should also be big, but if the value of property
    2 is big, the fitness value should be small.

    Giving a concrete example

    .. code-block:: python

        import stk

        normalizer = stk.Power((1, -1, 2))
        # Assuming that population holds molecule record instances
        # with the following fitness values
        # (1, 1, 1), (2, 2, 2), (3, 3, 3)
        # normalized will hold fitness values of
        # (1, 1, 1), (2, 0.5, 4), (3, 1/3, 9)
        normalized = tuple(normalizer.normalize(population))

    *Selectively Normalizing Fitness Values*

    Sometimes, you only want to normalize some members of a population,
    for example if some do not have an assigned fitness value,
    because the fitness calculation failed for whatever reason.
    You can use the `filter` parameter to exclude records from the
    normalization

    .. code-block:: python

        import stk

        normalizer = stk.Power(
            power=(1, 2, 3),
            # Only normalize values which are not None.
            filter=lambda population, record:
                record.get_fitness_value() is not None,
        )
        normalized = tuple(normalizer.normalize(population))

    """

    def __init__(
        self,
        power,
        filter=lambda population, record: True,
    ):
        """
        Initialize a :class:`.Power` instance.

        Parameters
        ----------
        power : :class:`float` or \
                :class:`tuple` of :class:`float`
            The power each fitness value is raised to. Can
            be a single number or multiple numbers, depending on the
            form of the fitness value.

        filter : :class:`callable`, optional
            Takes two parameters, first is a :class:`tuple`
            of :class:`.MoleculeRecord` instances,
            and the second is a :class:`.MoleculeRecord`. The
            :class:`callable` returns ``True`` or ``False``. Only
            molecules which return ``True`` will have fitness values
            normalized. By default, all molecules will have fitness
            values normalized.
            The instance passed to the `population` argument of
            :meth:`.normalize` is passed as the first argument, while
            the second argument will be passed every
            :class:`.MoleculeRecord` in it, one at a time.

        """

        self._power = power
        self._filter = filter

    def normalize(self, population):
        for record in population:
            if self._filter(population, record):
                yield record.with_fitness_value(
                    fitness_value=np.float_power(
                        record.get_fitness_value(),
                        self._power,
                    ),
                )
            else:
                yield record
