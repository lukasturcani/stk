"""
Fitness Normalizers
===================

#. :class:`.NullFitnessNormalizer`
#. :class:`.Power`
#. :class:`.Multiply`
#. :class:`.Sum`
#. :class:`.DivideByMean`
#. :class:`.ShiftUp`

Fitness normalizers are classes which are responsible for normalizing
the fitness values in a :class:`.Population`. They analyze the
:attr:`fitness` values across the entire population and update them.
Calling :meth:`~FitnessNormalizer.normalize` directly modifies the
:attr:`fitness` attribute of :class:`.Molecule` objects.

To see how :class:`FitnessNormalizer` can be used, look at the
documention of the classes which inherit it, for example
:class:`Power`, :class:`Sum` :class:`DivideByMean`. In addition,
multiple :class:`FitnessNormalizer` can be chained using
:meth:`.NormalizerSequence`.


.. _`adding fitness normalizers`:

Making New Fitness Normalizers
------------------------------

A new class inheriting :class:`FitnessNormalizer` must be made.
The class must define a :meth:`~FitnessNormalizer._normalize` method,
which takes one argument, which is a :class:`.Population` of
:class:`.Molecule`.


"""

import numpy as np
import logging


from ..base_calculators import Calculator


logger = logging.getLogger(__name__)


class FitnessNormalizer(Calculator):
    """
    Normalizes fitness values across a :class:`.Population`.

    """

    def __init__(self, filter=lambda mol: True, **kwargs):
        """
        Initialize a :class:`.FitnessNormalizer`.

        Parameters
        ----------
        filter : :class:`callable`, optional
            Takes a single parameter of type :class:`.Molecule` and
            returns ``True`` or ``False``. Only molecules which
            return ``True`` will have fitness values normalized. By
            default, all molecules will have fitness values normalized.

        """

        self._filter = filter
        super().__init__(filter=filter, **kwargs)

    def normalize(self, population, fitness_values):
        """
        Normalize the fitness values in `population`.

        Parameters
        ----------
        population : :class:`.iterable`
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

        normalized = dict(fitness_values)
        filtered = filter(self._filter, population)
        normalized.update(self._normalize(filtered, fitness_values))
        return normalized

    def _normalize(self, population, fitness_values):
        """
        Yield normalized `population` fitness values.

        Parameters
        ----------
        population : :class:`iterable`
            The filtered molecules.

        fitness_values : :class:`dict`
            Maps every molecule in `population` to its fitness value.

        Yields
        ------
        :class:`tuple`
            Holds a :class:`.Molecule` from `population` and its
            normalized value.

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method and needs to be implemented in a
            subclass.

        """

        raise NotImplementedError()


class NullFitnessNormalizer(FitnessNormalizer):
    """
    Does nothing.

    """

    def __init__(self, filter=lambda mol: True):
        """'
        Initialize a :class:`NullFitnessNormalizer`.

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
        return dict(fitness_values)


class Power(FitnessNormalizer):
    """
    Raises fitness values to some power.

    This works for cases where the :attr:`fitness` is single
    :class:`float` and where it is :class:`list` of :class:`float`.

    Examples
    --------
    Raising a :attr:`fitness` value by some power.

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
        power = stk.Power(2)

        # Normalize the fitness values.
        power.normalize(pop)

        # mol1.fitness is now 1.
        # mol2.fitness is now 4.
        # mol3.fitness is now 9

    Raising a :attr:`fitness` vector by some power.

    .. code-block:: python

        # Normally the fitness would be set by the get_fitness() method
        # of some a fitness calculator.
        mol1.fitness = [1, 2, 3]
        mol2.fitness = [4, 5, 6]
        mol3.fitness = [7, 8, 9]

        # Create the normalizer.
        power = stk.Power(2)

        # Normalize the fitness values.
        power.normalize(pop)

        # mol1.fitness is now [1, 4, 9].
        # mol2.fitness is now [16, 25, 36].
        # mol3.fitness is now [49, 64, 81]

    Raising a :attr:`get_fitness` vector by different powers.

    .. code-block:: python

        # Give the molecules some arbitrary fitness vectors.
        # Normally the fitness would be set by the get_fitness() method
        # of some a fitness calculator.
        mol1.fitness = [1, 2, 3]
        mol2.fitness = [4, 5, 6]
        mol3.fitness = [7, 8, 9]

        # Create the normalizer.
        power = stk.Power([1, 2, 3])

        # Normalize the fitness values.
        power.normalize(pop)

        # mol1.fitness is now [1, 4, 27].
        # mol2.fitness is now [4, 25, 216].
        # mol3.fitness is now [7, 64, 729]

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
        super().__init__(filter=filter)

    def _normalize(self, population, fitness_values):
        normalized = {}
        for mol in population:
            normalized[mol] = (
                np.float_power(fitness_values[mol], self.power)
            )
        return normalized


class Multiply(FitnessNormalizer):
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
        normalized = {}
        for mol in population:
            normalized[mol] = (
                np.multiply(fitness_values[mol], self._coefficient)
            )
        return normalized


class Sum(FitnessNormalizer):
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
        normalized = {}
        for mol in population:
            normalized[mol] = sum(fitness_values[mol])
        return normalized


class DivideByMean(FitnessNormalizer):
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

        normalized = {}
        for mol in population:
            normalized[mol] = np.divide(fitness_values[mol], mean)
        return normalized


class ShiftUp(FitnessNormalizer):
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

        # Get the minimum values of each element in the population.
        mins = np.min(fmat, axis=0)
        # Convert all the ones which are not to be shifted to 0 and
        # multiply the which are to be shifted by 1.01.
        is_array = isinstance(mins, np.ndarray)
        if not is_array:
            mins = np.array([mins])
        shift = np.zeros(len(mins))
        for i, min_ in enumerate(mins):
            if min_ <= 0:
                shift[i] = 1 - min_

        normalized = {}
        for mol in population:
            normalized[mol] = fitness_values[mol] + shift
        return normalized


class ReplaceFitness(FitnessNormalizer):
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
