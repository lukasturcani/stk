"""
Multiply
========

"""

import numpy as np
from functools import partial

from .fitness_normalizer import FitnessNormalizer


class Multiply(FitnessNormalizer):
    """
    Multiplies the fitness values by some coefficient.

    Examples
    --------
    *Multiplying Fitness Values by a set of Coefficients*

    In this example, assume that each fitness value consists of a
    :class:`tuple` of numbers, each representing a different property
    of the molecule, and each contributing to the final fitness value.
    The properties can be anything, such as  energy, number of atoms
    or diameter.

    If your final fitness value depends on the combination of these
    properties, you will probably first want to scale them with
    :class:`.DivideByMean`. Once this is done, you may want to
    multiply each property by some coefficient, which reflects its
    relative importance to the final fitness value. For example
    if you multiply the value of one property by ``1`` and another
    by ``2``, the second will contribute twice as much to the
    final fitness value, assuming that you get the final fitness value
    by using the :class:`.Sum` normalizer after :class:`.Multiply`.

    Giving a concrete example

    .. code-block:: python

        import stk

        normalizer = stk.Coefficient((1, 2, 3))
        # Assuming that population holds molecule record instances
        # with the following fitness values
        # (1, 1, 1), (2, 2, 2), (3, 3, 3)
        # normalized will hold fitness values of
        # (1, 2 ,3), (2, 4, 6), (3, 6, 9)
        normalized = tuple(normalizer.normalize(population))

    *Selectively Normalizing Fitness Values*

    Sometimes, you only want to normalize some members of a population,
    for example if some do not have an assigned fitness value,
    because the fitness calculation failed for whatever reason.
    You can use the `filter` parameter to exclude records from the
    normalization

    .. code-block:: python

        import stk

        normalizer = stk.Coefficient(
            coefficient=(1, 2, 3),
            # Only normalize values which are not None.
            filter=lambda population, record:
                record.get_fitness_value() is not None,
        )
        normalized = tuple(normalizer.normalize(population))



    """

    def __init__(
        self,
        coefficient,
        filter=lambda population, mol: True,
    ):
        """
        Initialize a :class:`Multiply` instance.

        Parameters
        ----------
        coefficient : :class:`float` or :class:`list` of :class:`float`
            The coefficients each :attr:`fitness` value by. Can be
            a single number or multiple numbers.

        filter : :class:`callable`, optional
            Takes a two parameters, first is a :class:`.EAPopulation`
            and the second is a :class:`.Molecule`, and
            returns ``True`` or ``False``. Only molecules which
            return ``True`` will have fitness values normalized. By
            default, all molecules will have fitness values normalized.
            The :class:`.EAPopulation` on which :meth:`normalize` is
            called is passed as the first argument while the second
            argument will be passed every :class:`.Molecule` in it.

        """

        self._coefficient = coefficient
        self._filter = filter

    def normalize(self, population):
        filtered = tuple(filter(
            partial(self._filter, population),
            population,
        ))
        for record in filtered:
            # Use np.multiply here so that both lists and arrays work
            # properly. If * is used [8]*3 will wrongly make [8, 8, 8].
            yield record.with_fitness_value(
                fitness_value=np.multiply(
                    self._coefficient,
                    record.get_fitness_value(),
                ),
            )
