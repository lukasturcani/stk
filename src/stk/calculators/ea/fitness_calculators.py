"""
Fitness Calculators
===================

#. :class:`.PropertyVector`
#. :class:`.RaisingFitnessCalculator`

Fitness calculators are classes which inherit
:class:`FitnessCalculator` and define a
:meth:`~FitnessCalculator.get_fitness` method. This method is
used to calculate the fitness of molecules. A
:class:`FitnessCalculator` will hold calculated fitness values in a
cache. The method will also create a :attr:`fitness` attribute on the
molecules it evaluates, which holds the fitness value. The values
calculated by :meth:`~FitnessCalculator.get_fitness` can be any Python
object, as long as the :attr:`fitness` value after
:class:`FitnessNormalizer.normalize` is
applied is set to a positive, non-zero :class:`float`. The calculated
fitness may also be ``None`` to indicate a failed calculation.

The :class:`FitnessCalculator` can be pickled if the calculated values
are to be saved.

For examples of how a :class:`FitnessCalculator` may be used, look
at the documentation of classes which inherit it, for example
:class:`PropertyVector`.

During the GA, fitness values are initially calculated by a
fitness calculator, which is an instance of a
:class:`FitnessCalculator`. After this, fitness normalization takes
place through an instance of a :class:`.FitnessNormalizer`. The
difference between fitness calculation and normalization is that
a fitness calculation always returns the same fitness value for a
given molecule and conformer, while fitness normalization updates
existing fitness values based on all the fitness values in a
population, for example by dividing the fitness value of all molecules
by the mean fitness across the population.


.. _`adding fitness calculators`:

Making New Fitness Calculators
------------------------------

A new class inheriting :class:`FitnessCalculator` must be created.
The class must define a :meth:`~FitnessCalculator.get_fitness` method,
which takes one parameter, `mol`, which takes a
:class:`.Molecule` object and is the molecule whose fitness is to be
calculated.

"""

from functools import wraps
import logging
import numpy as np

logger = logging.getLogger(__name__)


def _add_fitness_attribute_creation(get_fitness):
    """
    Make fitness functions add a :attr:`fitness` attribute.

    The attribute is added to the :class:`.Molecule` objects evaluated
    by the :meth:`~FitnessCalculator.get_fitness` method.

    Parameters
    ----------
    get_fitness : :class:`function`
        A fitness function, which is a
        :meth:`~FitnessCalculator.get_fitness` method of a
        :class:`FitnessCalculator`.

    Returns
    -------
    :class:`function`
        The decorated fitness function.

    """

    @wraps(get_fitness)
    def inner(self, mol):
        r = get_fitness(self, mol)
        mol.fitness = r
        return r

    return inner


def _add_cache_use(get_fitness):
    """
    Give fitness functions the option skip re-calculatations.

    Parameters
    ----------
    fitness : :class:`function`
        A fitness function, which is a
        :meth:`~FitnessCalculator.get_fitness` method of a
        :class:`FitnessCalculator`.

    Returns
    -------
    :class:`function`
        The decorated fitness function.

    """

    @wraps(get_fitness)
    def inner(self, mol):
        if self._use_cache and mol in self._cache:
            logger.info(
                f'Using cached fitness value for "{mol}".'
            )
            return self._cache[mol]
        else:
            r = get_fitness(self, mol)
            if self._use_cache:
                self._cache[mol] = r
            return r

    return inner


class FitnessCalculator:
    """
    Calculates and stores fitness values of molecules.

    A :class:`FitnessCalculator` will automatically add a
    :attr:`fitness` attribute to any :class:`.Molecule` objects
    it calculates a fitness value for. The attribute will hold the
    calculated fitness value

    """

    def __init__(self, use_cache):
        """
        Initialize a :class:`FitnessCalculator` instance.

        Parameters
        ----------
        use_cache : :class:`bool`
            If ``True`` then fitness values for molecules already held
            in the cache are not re-calculated but the value already
            stored is used.

        """

        self._use_cache = use_cache
        self._cache = {}

    def __init_subclass__(cls, **kwargs):
        cls.get_fitness = _add_cache_use(cls.get_fitness)
        cls.get_fitness = _add_fitness_attribute_creation(
            cls.get_fitness
        )
        return super().__init_subclass__(**kwargs)

    def set_cache_use(self, use_cache):
        """
        Set cache use on or off.

        Parameters
        ----------
        use_cache : :class:`bool`
            ``True`` if the cache is to be used.

        Returns
        -------
        None : :class:`NoneType`

        """

        self._use_cache = use_cache

    def is_caching(self):
        """
        ``True`` if the optimizer has caching turned on.

        Returns
        -------
        :class:`bool`
            ``True`` if the optimizer has caching turned on.

        """

        return self._use_cache

    def add_to_cache(self, mol, fitness):
        """
        Add the `fitness` of `mol` to the cache.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule whose `fitness` should be added to the cache.

        fitness : :class:`object`
            The fitness value to be added to the cache.

        Returns
        -------
        None : :class:`NoneType`

        """

        self._cache[mol] = fitness

    def get_fitness(self, mol):
        """
        Get the fitness value of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule whose fitness should be calculated.

        Returns
        -------
        :class:`object`
            The fitness value of `mol`. Can be
            ``None`` to indicate a failed calculation.

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method and needs to be implemented in a
            subclass.

        """

        raise NotImplementedError()


class RaisingFitnessCalculatorError(Exception):
    ...


class RaisingFitnessCalculator(FitnessCalculator):
    """
    Raises and calculates fitness at random.

    This fitness calculator is used for debugging to simulate
    cases where the fitness calculation raises an error.

    Examples
    --------
    .. code-block:: python

        import stk

        mol = stk.BuildingBlock('NCCNCCN', ['amine'])
        calc = stk.PropertyVector(lambda mol: len(mol.atoms))
        partial_raiser = RaisingFitnessCalculator(
            fitness_calculator=calc,
            fail_chance=0.75
        )
        # 75 % chance an error will be raised by calling get_fitness.
        partial_raiser.get_fitness(mol)

    """

    def __init__(
        self,
        fitness_calculator,
        fail_chance=0.5,
        use_cache=False
    ):
        """
        Initialize a :class:`RaisingFitnessCalculator` instance.

        Parameters
        ----------
        fitness_calculator : :class:`FitnessCalculator`
            When the fitness calculator does not fail, it uses this
            :class:`FitnessCalculator` to calculate fitness of
            molecules.

        fail_chance : :class:`float`, optional
            The probability that the fitness calculator will raise an
            error each time :meth:`fitness` is used.

        use_cache : :class:`bool`, optional
            If ``True`` then fitness values for molecules already held
            in the cache are not re-calculated but the value already
            stored is used.

        """

        self._fitness_calculator = fitness_calculator
        self._fail_chance = fail_chance
        super().__init__(use_cache=use_cache)

    def get_fitness(self, mol):
        """
        Get the fitness value of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule whose fitness should be calculated.

        Returns
        -------
        :class:`object`
            Whatever the :meth:`get_fitness` method of
            :attr:`fitness_calculator` returns.

        Raises
        ------
        :class:`RaisingFitnessCalculatorError`
            This error is raised at random.

        """

        if np.random.rand() < self._fail_chance:
            raise RaisingFitnessCalculatorError(
                    'Used RaisingFitnessCalculator'
            )
        return self._fitness_calculator.get_fitness(mol)


class PropertyVector(FitnessCalculator):
    """
    Calculates a set of properties of a molecule.

    This :class:`FitnessCalculator` applies a series of
    :class:`function` to a :class:`.Molecule` and appends each result
    to a :class:`list`. The :class:`list` forms the property vector of
    the molecule and it is returned as the fitness value of the
    molecule.

    If any of the :class:`function` returns ``None``, then instead of
    a :class:`list` the fitness value will be ``None``. In essence,
    :class:`PropertyVector` requires that all properties were
    successfully calculated.

    Examples
    --------
    Use on :class:`.BuildingBlock` objects.

    .. code-block:: python

        import stk

        # Create the molecules which are to have fitness values
        # evaluated.
        mol1 = stk.BuildingBlock('NCCN', ['amine'])
        mol2 = stk.BuildingBlock('NC[Si]CCCN', ['amine'])
        mol3 = stk.BuildingBlock('O=CCC(C=O)CCC=O', ['aldehyde'])

        # Create the functions which calculate the molecule properties.
        def atom_number(mol):
            return len(mol.atoms)

        def diameter(mol):
            return mol.maximum_diamater()

        def energy(mol):
            energy_calculator = stk.MMFFEnergy()
            return energy_calculator.get_energy(mol)

        # Create the fitness calculator.
        fitness_calculator = PropertyVector(
            atom_number,
            diameter,
            energy
        )

        # Calculate the fitness vector of mol1. It will be a list
        # holding the number of atoms, diameter and energy of mol1,
        # respectively.
        mol1_fitness = fitness_calculator.get_fitness(mol1)
        # The molecule will also have a fitness attribute holding the
        # result.
        if mol1.fitness == mol1_fitness:
            print('Fitness attribute added.')

        # Calculate the fitness vector of mol2. It will be a list
        # holding the number of atoms, diameter and energy of mol2,
        # respectively.
        mol2_fitness = fitness_calculator.get_fitness(mol2)
        # The molecule will also have a fitness attribute holding the
        # result.
        if mol2.fitness == mol2_fitness:
            print('Fitness attribute added.')

        # Calculate the fitness vector of mol3. It will be a list
        # holding the number of atoms, diameter and energy of mol3,
        # respectively.
        mol3_fitness = fitness_calculator.get_fitness(mol3)
        # The molecule will also have a fitness attribute holding the
        # result.
        if mol3.fitness == mol3_fitness:
            print('Fitness attribute added.')


    Use on :class:`.ConstructedMolecule` objects

    .. code-block:: python

        # First create molecules whose fitness value we wish to
        # caclculate.
        bb1 = stk.BuildingBlock('[Br]CC[Br]', ['bromine'])
        polymer1 = stk.ConstructedMolecule(
            building_blocks=[bb1],
            topology_graph=stk.polymer.Linear('A', [0], n=5)
        )

        bb2 = stk.BuildingBlock('[Br]CCNNCC[Br]', ['bromine'])
        polymer2 = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topology_graph=stk.polymer.Linear('AB', [0, 0], n=2)
        )

        # Create the functions which calculate the molecule properties.
        def atom_number(mol):
            return len(mol.atoms)

        def diameter(mol):
            return mol.maximum_diamater()

        def monomer_number(mol):
            return sum(mol.building_block_counter.values())

        # Create the fitness calculator.
        fitness_calculator = PropertyVector(
            atom_number,
            diameter,
            monomer_number
        )

        # Calculate the fitness vector of polymer1. It will be a list
        # holding the number of atoms, diameter and the number of
        # monomers in polymer1, espectively.
        polymer1_fitness = fitness_calculator.get_fitness(polymer1)
        # The molecule will also have a fitness attribute holding the
        # result.
        if polymer1.fitness == polymer1_fitness:
            print('Fitness attribute added.')

        # Calculate the fitness vector of polymer2. It will be a list
        # holding the number of atoms, diameter and the number of
        # monomers in polymer2, espectively.
        polymer2_fitness = fitness_calculator.get_fitness(polymer2)
        # The molecule will also have a fitness attribute holding the
        # result.
        if polymer2.fitness == polymer2_fitness:
            print('Fitness attribute added.')

    """

    def __init__(self, *property_fns, use_cache=False):
        """
        Initialize a :class:`CageFitness` instance.

        Parameters
        ----------
        *property_fns : :class:`tuple` of :class:`function`
            A group of :class:`function`, each of which is used to
            calculate a single property of the molecule. Each function
            must take one parameter, `mol`, which accepts
            a :class:`.Molecule` object. This is the molecule used to
            calculate the property.

        use_cache : :class:`bool`, optional
            If ``True`` then fitness values for molecules already held
            in the cache are not re-calculated but the value already
            stored is used.

        """

        self._property_fns = property_fns
        super().__init__(use_cache=use_cache)

    def get_fitness(self, mol):
        """
        Get the fitness of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule whose property vector should be calculated.

        Returns
        -------
        :class:`list`
            A :class:`list` of properties of the `mol`. Will be
            ``None`` if any of the properties fail.

        """

        property_vector = []
        for property_fn in self._property_fns:
            logger.info(
                f'Using {property_fn.__name__} on "{mol}".'
            )
            r = property_fn(mol)
            if r is None:
                logger.warning(
                    f'Using '
                    f'{property_fn.__name__} on "{mol}" failed.'
                )
            property_vector.append(r)
        if None in property_vector:
            return None
        else:
            return property_vector
