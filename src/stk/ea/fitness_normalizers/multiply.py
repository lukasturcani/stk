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
    Multiplying a fitness value by a coefficient.

    .. code-block:: python

        import stk

        # Create the normalizer.
        multiply = stk.Multiply(2)

        # Normalize the fitness values. Assume the fitness values are
        # {mol1: 1, mol2: 2, mol3: 3}.
        normalized = multiply.normalize(pop)

        # normalized is
        # {mol1: 2, mol2: 4, mol3: 6}.


    Multiplying a vector of fitness values by some coefficient

    .. code-block:: python

        multiply = stk.Multiply(2)

        # Normalize the fitness values. Assume the fitness values are
        # {mol1: [1, 2, 3], mol2: [4, 5, 6], mol3: [7, 8, 9]}.
        normalized = multiply.normalize(pop)

        # normalized is
        # {mol1: [2, 4, 6], mol2: [8, 10, 12], mol3: [14, 16, 18]}.


    Multiplying a vector of fitness values by different coefficients

    .. code-block:: python

        multiple = stk.Multiply([1, 2, 3])

        # Normalize the fitness values. Assume the fitness values are
        # {mol1: [1, 2, 3], mol2: [4, 5, 6], mol3: [7, 8, 9]}.
        normalized = multiply.normalize(pop)

        # normalized is
        # {mol1: [1, 4, 9], mol2: [4, 10, 18], mol3: [7, 16, 27]}.

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
