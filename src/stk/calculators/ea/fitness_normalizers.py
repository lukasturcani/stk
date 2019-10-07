"""
Fitness Normalizers
===================

#. :class:`.NullFitnessNormalizer`
#. :class:`.Power`
#. :class:`.Multiply`
#. :class:`.Sum`
#. :class:`.DivideByMean`
#. :class:`.ShiftUp`
#. :class:`.If`
#. :class:`.TryCatch`
#. :class:`.Sequence`
#. :class:`.Random`
#. :class:`.RaisingCalculator`

Fitness normalizers are classes which are responsible for normalizing
the fitness values in an :class:`.EAPopulation`. They analyze the
fitness values across the entire population and calculate new ones.

To see how :class:`FitnessNormalizer` can be used, look at the
documention of the classes which inherit it, for example
:class:`.Power`, :class:`.Sum` :class:`.DivideByMean`. In addition,
multiple :class:`.FitnessNormalizer` can be chained using
:meth:`.Sequence`.


.. _`adding fitness normalizers`:

Making New Fitness Normalizers
------------------------------

A new class inheriting :class:`FitnessNormalizer` must be made.
This is an abstract base class so its virtual methods need to be
implemented.

"""

import numpy as np
import logging


from ..base_calculators import Calculator


logger = logging.getLogger(__name__)


class FitnessNormalizer(Calculator):
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


class _FilteringNormalizer(FitnessNormalizer):
    """
    Implements some of the :class:`.FitnessNormalizer` interface.

    """

    def _normalize(self, population, fitness_values):
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

        normalized = population.get_fitness_values()
        filtered = filter(self._filter, population)
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


class NullFitnessNormalizer(FitnessNormalizer):
    """
    Does nothing.

    """

    def _normalize(self, population):
        return population.get_fitness_values()


class Power(_FilteringNormalizer, FitnessNormalizer):
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

    def __init__(self, power, filter=lambda mol: True):
        """
        Initialize a :class:`Power` instance.

        Parameters
        ----------
        power : :class:`float` or :class:`list` of :class:`float`
            The power to raise each :attr:`fitness` value to. Can be
            a single number or multiple numbers.

        filter : :class:`callable`, optional
            Takes a single parameter of type :class:`.Molecule` and
            returns ``True`` or ``False``. Only molecules which
            return ``True`` will have fitness values normalized. By
            default, all molecules will have fitness values normalized.

        """

        self._power = power
        self._filter = filter

    def _get_normalized_values(self, filtered, fitness_values):
        for mol in filtered:
            yield mol, np.float_power(fitness_values[mol], self._power)


class Multiply(_FilteringNormalizer, FitnessNormalizer):
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

    def __init__(self, coefficient, filter=lambda mol: True):
        """
        Initialize a :class:`Multiply` instance.

        Parameters
        ----------
        coefficient : :class:`float` or :class:`list` of :class:`float`
            The coefficients each :attr:`fitness` value by. Can be
            a single number or multiple numbers.

        filter : :class:`callable`, optional
            Takes a single parameter of type :class:`.Molecule` and
            returns ``True`` or ``False``. Only molecules which
            return ``True`` will have fitness values normalized. By
            default, all molecules will have fitness values normalized.

        """

        self._coefficient = coefficient
        self._filter = filter

    def _get_normalized_values(self, filtered, fitness_values):
        for mol in filtered:
            yield mol, self._coefficient*fitness_values[mol]


class Sum(_FilteringNormalizer, FitnessNormalizer):
    """
    Sums the values in a :class:`list`.

    Examples
    --------
    .. code-block:: python

        # Create the normalizer.
        sum_normalizer = stk.Sum()

        # Normalize the fitness values. Assume the fitness values are
        # {mol1: [1, 2, 3], mol2: [4, 5, 6], mol3: [7, 8, 9]}.
        normalized = sum_normalizer.normalize(pop)

        # normalized is
        # {mol1: 6, mol2: 15, mol3: 24}

    """

    def __init__(self, filter=lambda mol: True):
        """
        Initialize a :class:`.Sum` instance.

        Parameters
        ----------
        filter : :class:`callable`, optional
            Takes a single parameter of type :class:`.Molecule` and
            returns ``True`` or ``False``. Only molecules which
            return ``True`` will have fitness values normalized. By
            default, all molecules will have fitness values normalized.

        """

        self._filter = filter

    def _get_normalized_values(self, filtered, fitness_values):
        for mol in filtered:
            yield mol, sum(fitness_values[mol])


class DivideByMean(_FilteringNormalizer, FitnessNormalizer):
    """
    Divides fitness values by the population mean.

    While this function can be used if the fitness value of each
    :class:`.Molecule` in the :class:`.EAPopulation` is a single
    number, it is most useful when the fitness values is a
    :class:`list` of numbers. In this case, it is necessary to somehow
    combine the numbers so that a single fitness value is produced.
    For example, take a fitness value which is the vector holding the
    properties ``[energy, diameter, num_atoms]``. For a given molecule
    these numbers may be something like ``[200,000, 12, 140]``. If we
    were to sum these numbers, the energy term would dominate the final
    fitness value. In order to combine these numbers we can divide them
    by the population averages. For example, if the average energy
    of molecules in the population is ``300,000`` the average diameter
    is ``10`` and the average number of atoms is ``70`` then the
    fitness vector would be scaled to ``[0.5, 1.2, 2]``. These
    numbers are now of a similar magnitude and can be summed to give a
    reasonable value. After division , each value represents how
    much better than the population average each property value is.
    In essence we have removed the units from each parameter.

    Examples
    --------
    Scale fitness values

    .. code-block:: python

        import stk

        mean_scaler = stk.DivideByMean()

        # Normalize the fitness values.
        # Assume the fitness values are
        # {mol1: 1, mol2: 2, mol3: 3}
        normalized = mean_scaler.normalize(pop)

        # normalized is
        # {mol1: 0.5, mol2: 1, mol3: 1.5}


    Scale fitness vectors

    .. code-block:: python

        # Create the normalizer.
        # mean_scaler = DivideByMean()

        # Normalize the fitness values.
        # Assume the fitness values are
        # {mol1: [1, 10, 100], mol2: [2, 20, 100], mol3: [3, 30, 100]}.
        normalized = mean_scaler.normalize(pop)

        # normalized is
        # {
        #     mol1: [0.5, 0.5, 0.5],
        #     mol2: [1, 1, 1],
        #     mol3: [1.5, 1.5, 1.5]
        # }.

    """

    def __init__(self, filter=lambda mol: True):
        """
        Initialize a :class:`.DivideByMean` instance.

        Parameters
        ----------
        filter : :class:`callable`, optional
            Takes a single parameter of type :class:`.Molecule` and
            returns ``True`` or ``False``. Only molecules which
            return ``True`` will have fitness values normalized. By
            default, all molecules will have fitness values normalized.

        """

        self._filter = filter

    def _get_normalized_values(self, filtered, fitness_values):
        # filtered gets iterated through multiple times.
        filtered = list(filtered)
        mean = np.mean(
            a=[fitness_values[mol] for mol in filtered],
            axis=0,
        )
        logger.debug(f'Means used in DivideByMean: {mean}')

        for mol in filtered:
            yield mol, fitness_values[mol] / mean


class ShiftUp(_FilteringNormalizer, FitnessNormalizer):
    """
    Shifts negative fitness values to be positive.

    Assume you have a vector-valued fitness value, where each number
    represents a different property of the molecule

    .. code-block:: python

        {mol1: [1, -10, 1]}

    One way to convert the vector-valued fitness value into a
    scalar fitness value is by summing the elements, and the result in
    this case would be ``-8``. Clearly this doesn't work, because the
    resulting fitness value is not a positive number. To fix this,
    the ``-10`` should be shifted to a positive value.

    :class:`.ShiftUp` finds the minimum value of each element in the
    vector-valued fitness value across the entire population, and for
    element where this minimum value is less than ``0``, shifts up
    the element value for every molecule in the population, so that the
    minimum value in the entire population is ``1``.

    For example, take a population with the vector-valued fitness
    values

    .. code-block:: python

        fitness_values = {
            mol1: [1, -5, 5],
            mol2: [3, -10, 2],
            mol3: [2, 20, 1],
        }

    After normalization the fitness values will be.

    .. code-block:: python

        normalized  = {
            mol1: [1, 6, 5],
            mol2: [3, 1, 2],
            mol3: [2, 31, 1],
        }

    This :class:`.ShiftUp` also works when the fitness value is a
    single value.

    Examples
    --------
    .. code-block:: python

        # Create the normalizer.
        shifter = ShiftUp()

        # Normalize the fitness values. Assume the fitness values are
        # {mol1: [1, -2, 3], mol2: [4, 5, -6], mol3: [7, 8, 9]}.
        normalized = shifter.normalize(pop)

        # normalized is
        # {mol1: [1, 1, 10], mol2: [4, 8, 1], mol3: [7, 11, 16]}.

    """

    def __init__(self, filter=lambda mol: True):
        """
        Initialize a :class:`.ShiftUp` instance.

        Parameters
        ----------
        filter : :class:`callable`, optional
            Takes a single parameter of type :class:`.Molecule` and
            returns ``True`` or ``False``. Only molecules which
            return ``True`` will have fitness values normalized. By
            default, all molecules will have fitness values normalized.

        """

        self._filter = filter

    def _get_normalized_values(self, filtered, fitness_values):
        # filtered is iterated through multiple times.
        filtered = list(filtered)

        # Get all the fitness arrays in a matrix.
        fmat = np.array([fitness_values[mol] for mol in filtered])

        # Get the minimum value of each element across the population.
        # keepdims ensures that np.min returns a 1-D array, because
        # it will be True if fitness values are scalar and False if
        # they are array-valued.
        mins = np.min(fmat, axis=0, keepdims=len(fmat.shape) == 1)

        # Convert all elements in mins which are not to be shifted to 0
        # and make the shift equal to the minimum value + 1.
        shift = np.zeros(len(mins))
        for i, min_ in enumerate(mins):
            if min_ <= 0:
                shift[i] = 1 - min_

        for mol in filtered:
            yield mol, fitness_values[mol] + shift


class ReplaceFitness(_FitnessNormalizer, FitnessNormalizer):
    def __init__(self, replacement_fn, filter=lambda mol: True):
        """
        Initialize a :class:`.ReplaceFitness` instance.

        Parameters
        ----------
        replacement_fn : :class:`callable`
            Takes a single parameter, the :class:`.Population` which
            needs to be normalized, before it is filtered, and
            returns an :class:`object` which is used as the new
            fitness value for all molecules which pass the
            `filter`.

        filter : :class:`callable`, optional
            Takes a single parameter of type :class:`.Molecule` and
            returns ``True`` or ``False``. Only molecules which
            return ``True`` will have fitness values replaced. By
            default, all molecules will have fitness values replaced.

        """

        super().__init__(filter=filter)

    def normalize(self, population, fitness_values):
        """
        Normalize the fitness values in `population`.

        Parameters
        ----------
        population : :class:`.Population`
            The molecules which need to have their fitness values
            normalized.

        fitness_values : :class:`dict`
            Maps every molecule in `population` to its fitness value.

        Returns
        -------
        :class:`dict`
            Maps every molecule in `population` to its normalized
            fitness value.

        """

        replacement_value = self._replacement_fn(population)
        normalized = dict(fitness_values)
        for mol in filter(self._filter, population):
            normalized[mol] = replacement_value
        return normalized
