"""
Fitness Calculators
===================

#. :class:`.PropertyVector`
#. :class:`.RaisingFitnessCalculator`
#. :class:`.If`
#. :class:`.TryCatch`
#. :class:`.Random`
#. :class:`.RaisingCalculator`

Fitness calculators are classes which inherit
:class:`FitnessCalculator` and define a
:meth:`~FitnessCalculator.get_fitness` method. This method is
used to calculate the fitness of molecules. A
:class:`FitnessCalculator` will hold calculated fitness values in a
cache. The values calculated by :meth:`~FitnessCalculator.get_fitness`
can be any Python object, as long as the value after
:class:`FitnessNormalizer.normalize` is
applied is set to a positive, non-zero :class:`float`.
The :class:`FitnessCalculator` can be pickled if the calculated values
are to be saved.

For examples of how a :class:`FitnessCalculator` may be used, look
at the documentation of classes which inherit it, for example
:class:`PropertyVector`.

During the EA, fitness values are initially calculated by a
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
:class:`FitnessCalculator` is an abstract base class and its virtual
methods must be implemented.

"""

import logging

from ..base_calculators import MoleculeCalculator, _MoleculeCalculator


logger = logging.getLogger(__name__)


class FitnessCalculator(MoleculeCalculator):
    """
    Calculates fitness values of molecules.

    """

    def get_fitness(self, mol):
        """
        Return the fitness value of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule whose fitness value should be calculated.

        Returns
        -------
        :class:`object`
            The fitness value of `mol`.

        """

        return self._cache_result(self._get_fitness, mol)

    def _get_fitness(self, mol):
        """
        Return the fitness value of `mol`.

        Parameters
        ----------
        mol : :class:`.Molecule`
            The molecule whose fitness should be calculated.

        Returns
        -------
        :class:`object`
            The fitness value of `mol`.

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method and needs to be implemented in a
            subclass.

        """

        raise NotImplementedError()


class FitnessFunction(_MoleculeCalculator, FitnessCalculator):
    """
    Takes a function and uses it as a calculator.

    """

    def __init__(self, fitness_fn, use_cache=False):
        """
        Initialize a :class:`.FitnessFunction` instance.

        Parameters
        ----------
        fitness_fn : :class:`callable`
            Take a single parameter, the :class:`.Molecule` whose
            fitness needs to be calculated, and returns its
            fitness value.

        use_cache : :class:`bool`, optional
            If ``True`` a fitness calculation will not be performed on
            the same molecule twice, instead the previously returned
            value will be returned.

        """

        self._fitness_fn = fitness_fn
        super().__init__(use_cache=use_cache)

    def _get_fitness(self, mol):
        return self._fitness_fn(mol)


class PropertyVector(_MoleculeCalculator, FitnessCalculator):
    """
    Calculates a set of properties of a molecule.

    This :class:`FitnessCalculator` applies a series of
    :class:`function` to a :class:`.Molecule` and appends each result
    to a :class:`list`. The :class:`list` forms the property vector of
    the molecule and it is returned as the fitness value of the
    molecule.

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

        # Calculate the fitness vector of mol2. It will be a list
        # holding the number of atoms, diameter and energy of mol2,
        # respectively.
        mol2_fitness = fitness_calculator.get_fitness(mol2)

        # Calculate the fitness vector of mol3. It will be a list
        # holding the number of atoms, diameter and energy of mol3,
        # respectively.
        mol3_fitness = fitness_calculator.get_fitness(mol3)


    Use on :class:`.ConstructedMolecule` objects

    .. code-block:: python

        # First create molecules whose fitness value we wish to
        # calculate.
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
            return mol.maximum_diameter()

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
        # monomers in polymer1, respectively.
        polymer1_fitness = fitness_calculator.get_fitness(polymer1)

        # Calculate the fitness vector of polymer2. It will be a list
        # holding the number of atoms, diameter and the number of
        # monomers in polymer2, respectively.
        polymer2_fitness = fitness_calculator.get_fitness(polymer2)

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

    def _get_fitness(self, mol):
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
            property_vector.append(property_fn(mol))
        return property_vector
