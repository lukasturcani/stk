"""
Module for defining fitness functions.

How fitness values are calculated.
----------------------------------

The basic outline of the fitness calculation procedure for a single
molecule is as follows,

.. code-block:: python

    # 1. Start with a MacroMolecule instance.
    mol = Polymer([monomer1, monomer2], SomeTopology())

    # This new instance has had no fitness functions applied to it
    # and therefore its "fitness" and "unscaled_fitness" attributes
    # are empty.
    mol.fitness  # None
    mol.unscaled_fitness  # {}

    # 2. Apply a fitness function to the MacroMolecule.
    mol.unscaled_fitness['random_fitness'] = random_fitness(mol)

    # This adds the unscaled fitness value into the "unscaled_fitness"
    # attribute, but does not change the "fitness" attribute.
    mol.fitness  # None
    mol.unscaled_fitness  # {'random_fitness': 34}

    # 3.Copy the value from "unscaled_fitness" into "fitness".
    mol.fitness = deepcopy(mol.unscaled_fitness['random_fitness'])

    mol.fitness  # 34
    mol.unscaled_fitness  # {'random_fitness': 34}

    # 4. Apply any number of normalization functions to the fitness
    #    value. Normalization functions can be skipped altogether.

    # Normalization functions act on a population of MacroMolecules,
    # i.e. they normalize the fitness values across the population. So
    # first step is to put the molecule into a population.
    pop = GAPopulation(mol)

    # Apply the normalization function.
    norm_func(pop)

    mol.fitness  # 12
    mol.unscaled_fitness  # {'random_fitness': 34}

    # Multiple normalization functions may be applied.
    norm_func2(pop)

    mol.fitness  # 200
    mol.unscaled_fitness  # {'random_fitness': 34}

When running the GA, this process is much simplified,

.. code-block:: python

    # Start with a population of molecules, representing a GA
    # generation.
    pop = GAPopulation(mol1, mol2, mol3)

    mol1.fitness  # None
    mol1.unscaled_fitness  # {}

    mol2.fitness  # None
    mol2.unscaled_fitness  # {}

    mol3.fitness  # None
    mol3.unscaled_fitness  # {}

    # Calculate the unscaled fitness values of all molecules.
    pop.calculate_member_fitness()

    mol1.fitness  # None
    mol1.unscaled_fitness  # {'random_fitness': 12}

    mol2.fitness  # None
    mol2.unscaled_fitness  # {'random_fitness': 34}

    mol3.fitness  # None
    mol3.unscaled_fitness  # {'random_fitness': 22}

    # Apply the normalization functions.
    pop.normalize_fitness_values()

    mol1.fitness  # 3
    mol1.unscaled_fitness  # {'random_fitness': 12}

    mol2.fitness  # 1
    mol2.unscaled_fitness  # {'random_fitness': 34}

    mol3.fitness  # 2
    mol3.unscaled_fitness  # {'random_fitness': 22}

The method :meth:`.GAPopulation.normalize_fitness_values` automatically
copies the fitness value from :attr:`.MacroMolecule.unscaled_fitness`
into :attr:`.MacroMolecule.fitness` before applying any normalization
functions. A fitness function is only applied once per molecule,
because the result depends only on the molecule. However, normalization
functions are re-applied every generation. This is because for some
normalization functions, the returned value for a molecule may also
depend on the other molecules in the population. At each generation,
a fresh copy of the value from :attr:`.MacroMolecule.unscaled_fitness`
is made into :attr:`.MacroMolecule.fitness`. This means that the
normalization functions always start from the unscaled fitness value.

.. _`adding fitness functions`:

Extending stk: Adding fitness functions.
----------------------------------------

To add a new fitness function simply write it as a function in this
module. It will need to take the :class:`.MacroMolecule` instance as
its first argument and this argument should be called `macro_mol`. The
purpose of this is to help users identify which arguments are handled
automatically by ``stk`` and which they need to define in the input
file. The convention is that if the fitness function takes an argument
called `macro_mol`, they do not have to specify that argument in the
input file.

A fitness function must return the value which represents the fitness
of the molecule `macro_mol`. If a fitness function is meant to be
paired with a normalization function it can return any value or
object it likes. Just as long as the normalization functions know how
to deal with it and convert it to a number.

.. code-block:: python

    def random_fitness(macro_mol):
        return abs(np.random.normal(50, 20))

Obviously, this is a just toy example but all the key components of
fitness functions are there. More complex fitness functions can take
an arbitrary number arguments, and will likely use molecular properties
to calculate a fitness value. Here is a fitness function which rewards
molecules for having a desired size.

.. code-block:: python

    def size(macro_mol, desired_size):
        return abs(macro_mol.max_diameter() - desired_size)

A fitness function may be complex and may not fit neatly into a single
function. For example, the :func:`cage_target` needs to call
:func:`_generate_complexes` in order to sample various conformations
before outputting a fitness value. This is fine. Define helper
functions such as :func:`_generate_complexes` within this module but
make sure they are private. This means that names of helper functions
begin with a leading underscore.

.. _plotting-note:

A note on plotting.
-------------------

As mentioned before, some fitness functions may be complex and as a
result manipulate all sorts of data. Typically, in order to measure the
progress of a GA, the fitness values in the population are tracked
across generations. However, let's say that some hypothetical fitness
function also calculates the energies of molecules. It may be quite
interesting plot the evolution of energies across generations too, even
if the energy is not directly reflect in the final fitness value. If
this is the case, the fitness function may assign to the
:attr:`.MacroMolecule.progress_params` attribute of `macro_mol`,

.. code-block:: python

    def example_func(macro_mol):
        mol_energy = macro_mol.energy.some_energy_func()
        macro_mol.progress_params['example_func'] = [mol_energy]
        ...

Now a plot showing the change in ``mol_energy`` across generations will
be made too, along with the plot showing the changes in fitness.

What if two things are needed to be kept track of? Simple,

.. code-block:: python

    def example_func(macro_mol):
        mol_energy = macro_mol.energy.some_energy_func()
        mol_radius = macro_mol.max_diamater() / 2
        macro_mol.progress_params['example_func'] = [mol_energy,
                                                     mol_radius]
        ...

Great, now a separate progress plot ``mol_energy`` and ``mol_radius``
will be made.

How will the y-axes be labelled in each plot? The decorator
:func:`_param_labels` exists for this.

Let's create a basic outline of a some fitness function:

.. code-block:: python

    @_param_labels('Molecule Energy / J mol-1', 'Mol Radius / m-9')
    def some_fitness_fn(macro_mol, some_param):
        ...
        calculate_stuff()
        ...
        macro_mol.progress_params['some_fitness_fn'] = [mol_energy,
                                                        mol_radius]
        ...
        return fitness_value

If this function is used in the GA, a progress plot will be made for
each variable placed in :attr:`.MacroMolecule.progress_params` and
they will have their y-axes labelled ``'Molecule Energy / J mol-1'``
and ``'Molecule Radius / m-9'``, respectively.

"""

from functools import wraps
import logging

logger = logging.getLogger(__name__)


def _add_fitness_update(fitness):
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


def _add_fitness_caching(fitness):
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
        key = (mol, conformer)
        if self.use_cache and key in self.fitness_values:
            return self.fitness_values[(mol, conformer)]

        r = fitness(self, mol, conformer)
        self.fitness_values[key] = r
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
    use_cache : :bool`True`
        If ``True`` then fitness values for molecules and conformers
        already held in :attr:`fitness_values` are not re-calculated
        and the value stored is used.

    fitness_values : :class:`dict`
        Stores fitness values of molecules in the form:

        .. code-block:: python

            fitness_values = {
                (mol1, conf1): 12.2,
                (mol1, conf3): 124.31,
                (mol2, conf1): 0.2
            }

        where ``mol1`` and ``mol2`` are :class:`.Molecule` objects
        and ``conf1`` and ``conf3`` are :class:`int` which are the
        conformers used to calculate the fitness values.

    """

    def __init__(self, use_cache=True):
        self.use_cache = use_cache
        self.fitness_values = {}

    def __init_subclass__(cls, **kwargs):
        cls.fitness = _add_fitness_update(cls.fitness)
        cls.fitness = _add_fitness_caching(cls.fitness)

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
        :class:`float`
            The fitness value of a `conformer` of `mol`.

        """

        raise NotImplementedError()


# Provides labels for the progress plotter.
# @_param_labels('Cavity Difference',
#                'Window Difference',
#                'Asymmetry',
#                'Energy per Bond',
#                'Precursors Strain',
#                'Dihedral Strain')
class PropertyVector(FitnessCalculator):
    """
    Calculates the a set of properties of a molecule.

    This :class:`FitnessCalculator` applies a series of
    :class:`function`s to a :class:`.Molecule` and appends each result
    to a :class:`list`. The :class:`list` forms the property vector of
    the molecule and it is returned as the fitness value of the
    molecule.

    Attributes
    ----------
    property_fns : :class:`tuple` of :class:`function`
        A group of :class:`function`s, each of which is used to
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
        # its fitness_values attribute.
        print(fitness_calculator.fitness_values)


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
        # its fitness_values attribute.
        print(fitness_calculator.fitness_values)

    Use on :class:`.MacroMolecule` objects, :class:`.Cage`

    .. code-block:: python

        # First create molecules whose fitness value we wish to
        # caclculate.
        bb1 = StructUnit2.smiles_init('NCCN', ['amine'])
        bb2 = StructUnit3.smiles_init('O=CCCC(C=O)CC=O', ['aldehyde'])

        cage1 = Cage([bb1, bb2], FourPlusSix())
        cage2 = Cage([bb1, bb2], EightPlusTwelve())

        # Create the functions which calculate the molecule properties.
        def cavity_size(mol, conformer):
            return mol.cavity_size(conformer)

        def window_variance(mol, conformer):
            return mol.window_variance(conformer)

        # Create the fitness calculator.
        fitness_calculator = PropertyVector(cavity_size,
                                            window_variance)

        # Calculate the fitness vector of cage1. It will be a list
        # holding the cavity size and window variance, respectively.
        cage1_fitness = fitness_calculator.fitness(cage1)
        # The molecule will also have a fitness attribute holding the
        # result.
        if cage1.fitness == cage1_fitness:
            print('Fitness attribute added.')

        # Calculate the fitness vector of cage2. It will be a list
        # holding the cavity size and window variance, respectively.
        cage2_fitness = fitness_calculator.fitness(cage2)
        # The molecule will also have a fitness attribute holding the
        # result.
        if cage2.fitness == cage2_fitness:
            print('Fitness attribute added.')

        # The fitness calculate will have all the results saved in
        # its fitness_values attribute.
        print(fitness_calculator.fitness_values)


    """

    def __init__(self, *property_fns):
        """
        Initializes a :class:`CageFitness` instance.

        Parameters
        ----------
        *property_fns : :class:`tuple` of :class:`function`
            A group of :class:`function`s, each of which is used to
            calculate a single property of the molecule. Each function
            must take 2 arguments, `mol` and `conformer`. `mol` accepts
            a :class:`.Molecule` object and `conformer` accepts an
            :class:`int`. These are the molecule and the conformer id
            used to calculate the property.

        """

        self.property_fns = property_fns

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
            A :class:`list` of properties of the `mol`.

        """

        property_vector = []
        for property_fn in self.property_fns:
            logger.info(
                f'Using {property_fn.__name__} on "{mol.name}".'
            )
            property_vector.append(property_fn(mol, conformer))
        return property_vector
