"""
Module for defining fitness calculators.

Fitness calculators are classes which inherit
:class:`FitnessCalculator` and define a
:meth:`~FitnessCalculator.fitness` method. This method is used to
calculate the fitness of molecules. A :class:`FitnessCalculator` will
hold calculated fitness values in
:attr:`FitnessCalculator.cache`. The method will also
create a :attr:`fitness` attribute on the molecules it evaluates,
which holds the fitness value. The values calculated by
:meth:`~FitnessCalculator.fitness` can be any Python object, as long as
the :attr:`fitness` value after :class:`FitnessNormalizer.normalize` is
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

Extending stk: Making new fitness calculators.
----------------------------------------------

A new class inheriting :class:`FitnessCalculator` must be created.
The class must define a :meth:`~FitnessCalculator.fitness` method,
which takes two arguments. The first is `mol` which takes a
:class:`.Molecule` object and is the molecule whose fitness is to be
calculated. The second is `conformer` and it is an :class:`int` holding
the conformer id of the conformer used for calculating the fitness.
`conformer` is should be an optional argument, defaulting to ``-1``.


"""

from functools import wraps
import logging
import numpy as np

logger = logging.getLogger(__name__)


def _add_fitness_attribute_creation(fitness):
    """
    Makes fitness functions add a :attr:`fitness` attribute.

    The attribute is added to the :class:`.Molecule` objects evaluated
    by the :meth:`~FitnessCalculator.fitness` method.

    Parameters
    ----------
    fitness : :class:`function`
        A fitness function, which is a
        :meth:`~FitnessCalculator.fitness` method of a
        :class:`FitnessCalculator`.

    Returns
    -------
    :class:`function`
        The decorated fitness function.

    """

    @wraps(fitness)
    def inner(self, mol, conformer=-1):
        r = fitness(self, mol, conformer)
        mol.fitness = r
        return r

    return inner


def _add_cache_use(fitness):
    """
    Gives fitness functions the option skip re-calculatations.

    Parameters
    ----------
    fitness : :class:`function`
        A fitness function, which is a
        :meth:`~FitnessCalculator.fitness` method of a
        :class:`FitnessCalculator`.

    Returns
    -------
    :class:`function`
        The decorated fitness function.

    """

    @wraps(fitness)
    def inner(self, mol, conformer=-1):
        key = (mol.key, conformer)
        if self.use_cache and key in self.cache:
            logger.info(
                'Using cached fitness value '
                f'for "{mol.name}" conformer {conformer}.'
            )
            return self.cache[key]
        else:
            r = fitness(self, mol, conformer)
            if self.use_cache:
                self.cache[key] = r
            return r

    return inner


class FitnessCalculator:
    """
    Calculates and stores fitness values of molecules.

    A :class:`FitnessCalculator` will automatically add a
    :attr:`fitness` attribute to any :class:`.Molecule` objects
    it calculates a fitness value for. The attribute will hold the
    calculated fitness value

    Attributes
    ----------
    use_cache : :class:`bool`
        If ``True`` then fitness values for molecules and conformers
        already held in :attr:`cache` are not re-calculated
        and the value already stored is used.

    cache : :class:`dict`
        Stores fitness values of molecules in the form:

        .. code-block:: python

            cache = {
                (mol1, conf1): 12.2,
                (mol1, conf3): 124.31,
                (mol2, conf1): 0.2
            }

        where ``mol1`` and ``mol2`` are :class:`.Molecule` objects
        and ``conf1`` and ``conf3`` are :class:`int` which are the
        conformers used to calculate the fitness values.

    """

    def __init__(self, use_cache):
        """
        Initializes a :class:`FitnessCalculator` instance.

        Parameters
        ----------
        use_cache : :class:`bool`
            If ``True`` then fitness values for molecules and
            conformers already held in :attr:`cache` are not
            re-calculated and the value already stored is used.

        """

        self.use_cache = use_cache
        self.cache = {}

    def __init_subclass__(cls, **kwargs):
        cls.fitness = _add_cache_use(cls.fitness)
        cls.fitness = _add_fitness_attribute_creation(cls.fitness)
        return super().__init_subclass__(**kwargs)

    def fitness(self, mol, conformer=-1):
        """
        Calculates the fitness value of a molecule.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule whose fitness should be calculated.

        conformer : :class:`int`, optional
            The conformer of `mol` to use.

        Returns
        -------
        :class:`object`
            The fitness value of a `conformer` of `mol`. Can be
            ``None`` to indicate a failed calculation.

        """

        raise NotImplementedError()


class RaisingFitnessCalculatorError(Exception):
    ...


class RaisingFitnessCalculator(FitnessCalculator):
    """
    Raises and calculates fitness at random.

    This fitness calculator is used for debugging to simulate
    cases where the fitness calculation raises an error.

    Attributes
    ----------
    fitness_calculator : :class:`FitnessCalculator`
        When the fitness calculator does not fail, it uses this
        :class:`FitnessCalculator` to calculate fitness of molecules.

    fail_chance : :class:`float`
        The probability that the fitness calculator will raise an error
        each time :meth:`fitness` is used.

    Examples
    --------
    .. code-block:: python

        mol = StructUnit.smiles_init('NCCNCCN', ['amine'])
        calc = PropertyVector(lambda mol, conformer: mol.GetNumAtoms())
        partial_raiser = RaisingFitnessCalculator(calc,
                                                  fail_chance=0.75)
        # 75 % chance an error will be raised by calling fitness.
        partial_raiser.fitness(mol)

    """

    def __init__(self,
                 fitness_calculator,
                 fail_chance=0.5,
                 use_cache=False):
        """
        Initializes a :class:`RaisingFitnessCalculator` instance.

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
            If ``True`` then fitness values for molecules and
            conformers already held in :attr:`cache` are not
            re-calculated and the value already stored is used.

        """

        self.fitness_calculator = fitness_calculator
        self.fail_chance = fail_chance
        super().__init__(use_cache=use_cache)

    def fitness(self, mol, conformer=-1):
        """
        Calculates the fitness of `conformer` of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule whose fitness should be calculated.

        conformer : :class:`int`, optional
            The conformer of `mol` to use.

        Returns
        -------
        :class:`object`
            Whatever the :meth:`fitness` method of
            :attr:`fitness_calculator` returns.

        Raises
        ------
        :class:`RaisingFitnessCalculatorError`
            This error is raised at random.

        """

        if np.random.rand() < self.fail_chance:
            raise RaisingFitnessCalculatorError(
                    'Used RaisingFitnessCalculator'
            )
        return self.fitness_calculator.fitness(mol, conformer)


class PropertyVector(FitnessCalculator):
    """
    Calculates the a set of properties of a molecule.

    This :class:`FitnessCalculator` applies a series of
    :class:`function` to a :class:`.Molecule` and appends each result
    to a :class:`list`. The :class:`list` forms the property vector of
    the molecule and it is returned as the fitness value of the
    molecule.

    If any of the :class:`function` returns ``None``, then instead of
    a :class:`list` the fitness value will be ``None``. In essence,
    :class:`PropertyVector` requires that all properties were
    successfully calculated.

    Attributes
    ----------
    property_fns : :class:`tuple` of :class:`function`
        A group of :class:`function`, each of which is used to
        calculate a single property of the molecule. Each function must
        take 2 arguments, `mol` and `conformer`. `mol` accepts a
        :class:`.Molecule` object and `conformer` accepts an
        :class:`int`. These are the molecule and the conformer id used
        to calculate the property.


    Examples
    --------
    Use on :class:`.StructUnit` objects.

    .. code-block:: python

        # Create the molecules which are to have fitness values
        # evaluated.
        mol1 = StructUnit.smiles_init('NCCN', ['amine'])
        mol2 = StructUnit2.smiles_init('NC[Si]CCCN', ['amine'])
        mol3 = StructUnit3.smiles_init('O=CCC(C=O)CCC=O', ['aldehyde'])

        # Create the functions which calculate the molecule properties.
        def atom_number(mol, conformer):
            return mol.mol.GetNumAtoms()

        def diameter(mol, conformer):
            return mol.max_diamater(conformer)

        def energy(mol, conformer):
            energy_calculator = MMFFEnergy()
            return energy_calculator.energy(mol, conformer)

        # Create the fitness calculator.
        fitness_calculator = PropertyVector(atom_number,
                                            diameter,
                                            energy)

        # Calculate the fitness vector of mol1. It will be a list
        # holding the number of atoms, diameter and energy of mol1,
        # respectively.
        mol1_fitness = fitness_calculator.fitness(mol1)
        # The molecule will also have a fitness attribute holding the
        # result.
        if mol1.fitness == mol1_fitness:
            print('Fitness attribute added.')

        # Calculate the fitness vector of mol2. It will be a list
        # holding the number of atoms, diameter and energy of mol2,
        # respectively.
        mol2_fitness = fitness_calculator.fitness(mol2)
        # The molecule will also have a fitness attribute holding the
        # result.
        if mol2.fitness == mol2_fitness:
            print('Fitness attribute added.')

        # Calculate the fitness vector of mol3. It will be a list
        # holding the number of atoms, diameter and energy of mol3,
        # respectively.
        mol3_fitness = fitness_calculator.fitness(mol3)
        # The molecule will also have a fitness attribute holding the
        # result.
        if mol3.fitness == mol3_fitness:
            print('Fitness attribute added.')

        # The fitness calculate will have all the results saved in
        # its cache attribute.
        print(fitness_calculator.cache)


    Use on :class:`.MacroMolecule` objects, :class:`.Polymer`

    .. code-block:: python

        # First create molecules whose fitness value we wish to
        # caclculate.
        bb1 = StructUnit2.smiles_init('[Br]CC[Br]', ['bromine'])
        polymer1 = Polymer([bb1], Linear('A', [0], n=5))

        bb2 = StructUnit2.smiles_init('[Br]CCNNCC[Br]', ['bromine'])
        polymer2 = Polymer([bb1, bb2], Linear('AB', [0, 0], n=2))

        # Create the functions which calculate the molecule properties.
        def atom_number(mol, conformer):
            return mol.mol.GetNumAtoms()

        def diameter(mol, conformer):
            return mol.max_diamater(conformer)

        def monomer_number(mol, conformer):
            return mol.topology.n * len(mol.topology.repeating_unit)

        # Create the fitness calculator.
        fitness_calculator = PropertyVector(atom_number,
                                            diameter,
                                            monomer_number)

        # Calculate the fitness vector of polymer1. It will be a list
        # holding the number of atoms, diameter and the number of
        # monomers in polymer1, espectively.
        polymer1_fitness = fitness_calculator.fitness(polymer1)
        # The molecule will also have a fitness attribute holding the
        # result.
        if polymer1.fitness == polymer1_fitness:
            print('Fitness attribute added.')

        # Calculate the fitness vector of polymer2. It will be a list
        # holding the number of atoms, diameter and the number of
        # monomers in polymer2, espectively.
        polymer2_fitness = fitness_calculator.fitness(polymer2)
        # The molecule will also have a fitness attribute holding the
        # result.
        if polymer2.fitness == polymer2_fitness:
            print('Fitness attribute added.')

        # The fitness calculate will have all the results saved in
        # its cache attribute.
        print(fitness_calculator.cache)

    Use on :class:`.MacroMolecule` objects, :class:`.Cage`

    .. code-block:: python

        # First create molecules whose fitness value we wish to
        # caclculate.
        bb1 = StructUnit2.smiles_init('NCCN', ['amine'])
        bb2 = StructUnit3.smiles_init('O=CCCC(C=O)CC=O', ['aldehyde'])

        cage1 = Cage([bb1, bb2], FourPlusSix())
        cage2 = Cage([bb1, bb2], EightPlusTwelve())

        # Create the functions which calculate the molecule properties.
        def cavity_difference(mol, conformer):
            # Calculate the difference from a cavity size of 5.
            return abs(mol.cavity_size(conformer)-5)

        def window_variance(mol, conformer):
            return mol.window_variance(conformer)

        # Create the fitness calculator.
        fitness_calculator = PropertyVector(cavity_difference,
                                            window_variance)

        # Calculate the fitness vector of cage1. It will be a list
        # holding how different the cavity size is from 5 Angstrom and
        # window variance, respectively.
        cage1_fitness = fitness_calculator.fitness(cage1)
        # The molecule will also have a fitness attribute holding the
        # result.
        if cage1.fitness == cage1_fitness:
            print('Fitness attribute added.')

        # Calculate the fitness vector of cage2. It will be a list
        # holding how different the cavity size is from 5 Angstrom and
        # window variance, respectively.
        cage2_fitness = fitness_calculator.fitness(cage2)
        # The molecule will also have a fitness attribute holding the
        # result.
        if cage2.fitness == cage2_fitness:
            print('Fitness attribute added.')

        # The fitness calculate will have all the results saved in
        # its cache attribute.
        print(fitness_calculator.cache)

    """

    def __init__(self, *property_fns, use_cache=False):
        """
        Initializes a :class:`CageFitness` instance.

        Parameters
        ----------
        *property_fns : :class:`tuple` of :class:`function`
            A group of :class:`function`, each of which is used to
            calculate a single property of the molecule. Each function
            must take 2 arguments, `mol` and `conformer`. `mol` accepts
            a :class:`.Molecule` object and `conformer` accepts an
            :class:`int`. These are the molecule and the conformer id
            used to calculate the property.

        use_cache : :class:`bool`, optional
            If ``True`` then fitness values for molecules and
            conformers already held in :attr:`cache` are not
            re-calculated and the value stored is used.

        """

        self.property_fns = property_fns
        super().__init__(use_cache=use_cache)

    def fitness(self, mol, conformer=-1):
        """
        Returns the property vector of a molecule.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule whose property vector should be calculated.

        conformer : :class:`int`, optional
            The conformer of `mol` to use.

        Returns
        -------
        :class:`list`
            A :class:`list` of properties of the `mol`. Will be
            ``None`` if any of the properties fail.

        """

        property_vector = []
        for property_fn in self.property_fns:
            logger.info(
                f'Using {property_fn.__name__} on "{mol.name}".'
            )
            r = property_fn(mol, conformer)
            if r is None:
                logger.warning(
                    f'Using '
                    f'{property_fn.__name__} on "{mol.name}" failed.'
                )
            property_vector.append(r)
        if None in property_vector:
            return None
        else:
            return property_vector
