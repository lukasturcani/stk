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
            yield mol, np.float_power(fitness_values[mol], self.power)


class Multiply(_FilteringNormalizer, FitnessNormalizer):
    """
    Multiplies the fitness value by some coefficient.

    Examples
    --------
    Multiplying a :attr:`fitness` value by a coefficient.

    .. code-block:: python

        import stk

        # Create the molecules and give them some arbitrary fitness.
        # Normally the fitness would be set by the get_fitness() method
        # of some a fitness calculator.

        mol1 = stk.BuildingBlock('NCCCN')
        mol2 = stk.BuildingBlock('[Br]CCC[Br]', ['bromine'])
        mol3 = stk.ConstructedMolecule(
            building_blocks=[mol2],
            topology_graph=stk.polymer.Linear('A', [0], 5)
        )

        mol1.fitness = 1
        mol2.fitness = 2
        mol3.fitness = 3

        # Place the molecules in a Population.
        pop = stk.Population(mol1, mol2, mol3)

        # Create the normalizer.
        multiply = stk.Multiply(2)

        # Normalize the fitness values.
        multiply.normalize(pop)

        # mol1.fitness is now 2.
        # mol2.fitness is now 4.
        # mol3.fitness is now 6

    Multiplying a :attr:`fitness` vector by some coefficient.

    .. code-block:: python

        # Give the molecules some arbitrary fitness vectors.
        # Normally the fitness would be set by the get_fitness() method
        # of some a fitness calculator.
        mol1.fitness = [1, 2, 3]
        mol2.fitness = [4, 5, 6]
        mol3.fitness = [7, 8, 9]

        # Place the molecules in a Population.
        pop = stk.Population(mol1, mol2, mol3)

        # Create the normalizer.
        multiply = stk.Multiply(2)

        # Normalize the fitness values.
        multiply.normalize(pop)

        # mol1.fitness is now [2, 4, 6].
        # mol2.fitness is now [8, 10, 12].
        # mol3.fitness is now [14, 16, 18]

    Multiplying a :attr:`fitness` vector by different coefficients.

    .. code-block:: python

        # Give the molecules some arbitrary fitness vectors.
        # Normally the fitness would be set by the get_fitness() method
        # of some a fitness calculator.
        mol1.fitness = [1, 2, 3]
        mol2.fitness = [4, 5, 6]
        mol3.fitness = [7, 8, 9]

        # Create the normalizer.
        multiply = stk.Multiply([1, 2, 3])

        # Normalize the fitness values.
        multiply.normalize(pop)

        # mol1.fitness is now [1, 4, 9].
        # mol2.fitness is now [4, 10, 18].
        # mol3.fitness is now [7, 16, 27]

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
        super().__init__(filter=filter)

    def _normalize(self, population, fitness_values):
        for mol in population:
            yield mol, self._coefficient*fitness_values[mol]


class Sum(_FitnessNormalizer, FitnessNormalizer):
    """
    Sums the values in a :class:`list`.

    Examples
    --------
    .. code-block:: python

        mol1 = stk.BuildingBlock('NCCCN')
        mol2 = stk.BuildingBlock('[Br]CCC[Br]', ['bromine'])
        mol3 = stk.ConstructedMolecule(
            building_blocks=[mol2],
            topology_graph=stk.polymer.Linear('A', [0], 5)
        )

        # Create the molecules and give them some arbitrary fitness
        # vectors.
        # Normally the fitness would be set by the fitness() method of
        # some a fitness calculator.
        mol1.fitness = [1, 2, 3]
        mol2.fitness = [4, 5, 6]
        mol3.fitness = [7, 8, 9]

        # Place the molecules in a Population.
        pop = stk.Population(mol1, mol2, mol3)

        # Create the normalizer.
        sum_normalizer = stk.Sum()

        # Normalize the fitness values.
        sum_normalizer.normalize(pop)

        # mol1.fitness is now 6.
        # mol2.fitness is now 15.
        # mol3.fitness is now 24.

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

        super().__init__(filter=filter)

    def _normalize(self, population, fitness_values):
        for mol in population:
            yield mol, sum(fitness_values[mol])


class DivideByMean(_FitnessNormalizer, FitnessNormalizer):
    """
    Divides fitness values by the population mean.

    While this function can be used if the :attr:`fitness` attribute
    of each :class:`.Molecule` in the :class:`.Population` is a single
    number it is most useful when :attr:`fitness` is a :class:`list`
    of numbers. In this case, it is necessary to somehow combine the
    numbers so that a single :attr:`fitness` value is produced.
    For example, take a :attr:`fitness` vector holding the properties
    ``[energy, diameter, num_atoms]``. For a given molecule these
    numbers may be something like ``[200,000, 12, 140]``. If we were
    to sum these numbers, the energy term would dominate the final
    fitness value. In order to combine these numbers we can divide them
    by the population averages. For example, if the average energy
    of molecules in the population is ``300,000`` the average diameter
    is ``10`` and the average number of atoms is ``70`` then the
    fitness vector would be scaled to ``[0.5, 1.2, 2]``. These
    numbers are now of a similar magnitude and can be some to give a
    reasonable value. After scaling each parameter represents how
    much better than the population average each property value is.
    In essence we have removed the units from each parameter.

    Examples
    --------
    Scale fitness values.

    .. code-block:: python

        import stk

        # Create the molecules and give them some arbitrary fitness
        # vectors.
        # Normally the fitness would be set by the fitness() method of
        # some a fitness calculator.
        mol1 = stk.BuildingBlock('NCCCN')
        mol2 = stk.BuildingBlock('[Br]CCC[Br]', ['bromine'])
        mol3 = stk.ConstructedMolecule(
            building_blocks=[mol2],
            topology_graph=stk.polymer.Linear('A', [0], 5)
        )

        mol1.fitness = 1
        mol2.fitness = 2
        mol3.fitness = 3

        # Place the molecules in a Population.
        pop = stk.Population(mol1, mol2, mol3)

        # Create the normalizer.
        mean_scaler = stk.DivideByMean()

        # Normalize the fitness values.
        mean_scaler.normalize(pop)

        # mol1.fitness is now 0.5.
        # mol2.fitness is now 1.
        # mol3.fitness is now 1.5.

    Scale fitness vectors.

    .. code-block:: python

        # Give the molecules some arbitrary fitness vectors.
        # Normally the fitness would be set by the get_fitness() method
        # of some a fitness calculator.
        mol1.fitness = [1, 10, 100]
        mol2.fitness = [2, 20, 200]
        mol3.fitness = [3, 30, 300]

        # Create the normalizer.
        mean_scaler = DivideByMean()

        # Normalize the fitness values.
        mean_scaler.normalize(pop)

        # mol1.fitness is now [0.5, 0.5, 0.5].
        # mol2.fitness is now [1, 1, 1].
        # mol3.fitness is now [1.5, 1.5, 1.5].

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

        super().__init__(filter=filter)

    def _normalize(self, population, fitness_values):
        mean = np.mean(
            a=[fitness_values[mol] for mol in population],
            axis=0,
        )
        logger.debug(f'Means used in DivideByMean: {mean}')

        for mol in population:
            yield mol, fitness_values[mol] / mean


class ShiftUp(_FitnessNormalizer, FitnessNormalizer):
    """
    Shifts negative values to be positive.

    Assume you have a fitness vector, where each number represents
    a different property of the molecule

    .. code-block:: python

        mol.fitness = [1, -10, 1]

    One way to convert the fitness array into a fitness value is
    by summing the elements, and the result in this case would be
    ``-8``. Clearly this doesn't work, because the resulting fitness
    value is not a positive number. To fix this, the ``-10`` should be
    shifted to a positive value.

    This :class:`FitnessNormalizer` finds the minimum value of each
    property across the entire population, and for properties where
    this minimum value is less than ``0``, shifts up the property value
    for every molecule in the population, so that the minimum value is
    ``1``.

    For example, take a population with the fitness vectors

    .. code-block:: python

        mol1.fitness = [1, -5, 5]
        mol2.fitness = [3, -10, 2]
        mol3.fitness = [2, 20, 1]

    After normalization fitness vectors will be.

    .. code-block:: python

        mol1.fitness  # [1, 6, 5]
        mol2.fitness  # [3, 1, 2]
        mol3.fitness  # [2, 31, 1]

    This :class:`FitnessNormalizer` also works when the :attr:`fitness`
    is a single value.

    Examples
    --------
    .. code-block:: python

        # Create the molecules and give them some arbitrary fitness
        # vectors.
        # Normally the fitness would be set by the fitness() method of
        # some a fitness calculator.
        mol1 = stk.BuildingBlock('NCCCN')
        mol2 = stk.BuildingBlock('[Br]CCC[Br]', ['bromine'])
        mol3 = stk.ConstructedMolecule(
            building_blocks=[mol2],
            topology_graph=stk.polymer.Linear('A', [0], 5)
        )

        mol1.fitness = [1, -2, 3]
        mol2.fitness = [4, 5, -6]
        mol3.fitness = [7, 8, 9]

        # Place the molecules in a Population.
        pop = Population(mol1, mol2, mol3)

        # Create the normalizer.
        shifter = ShiftUp([1, 2, 3])

        # Normalize the fitness values.
        shifter.normalize(pop)

        # mol1.fitness is now [1, 1, 10].
        # mol2.fitness is now [4, 8, 1].
        # mol3.fitness is now [7, 11, 16]

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

        super().__init__(filter=filter)

    def _normalize(self, population, fitness_values):
        # Get all the fitness arrays in a matrix.
        fmat = np.array([fitness_values[mol] for mol in population])

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

        for mol in population:
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
