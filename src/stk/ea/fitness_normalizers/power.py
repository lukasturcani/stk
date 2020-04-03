"""
Power
=====

"""

import numpy as np
from functools import partial

from .fitness_normalizer import FitnessNormalizer


class Power(FitnessNormalizer):
    """
    Raises fitness values to some power.

    This works for cases where the fitness value is single
    :class:`float` and where it is :class:`list` of :class:`float`.

    Examples
    --------
    Raising a fitness value by some power

    .. code-block:: python

        import stk

        pop = stk.Population(...)
        # Assume this returns {mol1: 1, mol2: 2, mol3: 3}.
        pop.get_fitness_values()

        # Create the normalizer.
        power = stk.Power(2)

        # Normalize the fitness values.
        normalized = power.normalize(pop)

        # normalized is {mol1: 1, mol2: 4, mol3: 9}.


    Raising vector valued fitness values by some power

    .. code-block:: python

        # Create the normalizer.
        power = stk.Power(2)

        # Normalize the fitness values. Assume the fitness values are
        # {mol1: [1, 2, 3], mol2: [4, 5, 6], mol3: [7, 8, 9]}.
        normalized = power.normalize(pop)

        # normalized is
        # {mol1: [1, 4, 9], mol2: [16, 25, 36], mol3: [49, 64, 81]}.


    Raising vector valued fitness values by different powers

    .. code-block:: python

        # Create the normalizer.
        power = stk.Power([1, 2, 3])

        # Normalize the fitness values. Assume the fitness values are
        # {mol1: [1, 2, 3], mol2: [4, 5, 6], mol3: [7, 8, 9]}.

        # Normalize the fitness values.
        normalized = power.normalize(pop)

        # normalized is
        # {mol1: [1, 4, 27], mol2: [4, 25, 216], mol3: [7, 64, 729]}.

    """

    def __init__(self, power, filter=lambda population, record: True):
        """
        Initialize a :class:`.Power` instance.

        Parameters
        ----------
        power : :class:`float` or \
                :class:`tuple` of :class:`float`
            The power to raise each fitness value to. Can be
            a single number or multiple numbers, depending on the
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
        filtered = tuple(filter(
            partial(self._filter, population),
            population,
        ))
        for record in filtered:
            yield record.with_fitness_value(
                fitness_value=np.float_power(
                    record.get_fitness_value(),
                    self._power,
                ),
            )
