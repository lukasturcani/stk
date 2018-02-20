"""
Defines normalization functions via the :class:`.Normalization` class.

Normalization functions are functions that recalculate the fitness
values of members of a population. The difference between fitness
and normalization functions is that fitness functions only use the
:class:`.MacroMolecule` to calculate fitness. Normalzation functions
have access to all :class:`.MacroMolecules` in a population. As a
result, they can scale the fitness values across the entire population.
This is useful if you want to ensure a spread of fitness values in the
population.

At each generation, a number of normalization functions can be applied
in sequence.

.. _`adding normalization functions`:

Extending mtk: Adding normalization functions.
----------------------------------------------

If a new normalization function is to be added to ``mtk`` it should be
added as a method in :class:`.Normalization`. The only requirements are
that the first argument is `population` (excluding `self`).

The naming requirement exists to help users identify which arguments
are handled automatically by the GA and which they need to define in
the input file. The convention is that if the normalization function
takes an argument called  `population`, it does not have to be
specified in the input file.

Normalization functions should not interact with
:attr:`.MacroMolecule.unscaled_fitness` in any way. They should only
modify the value in :attr:`.MacroMolecule.fitness`. Before the first
normalization function is applied by a :class:`Normalization` instance,
the value in :attr:`~.MacroMolecule.unscaled_fitness` is copied into
:attr:`~.MacroMolecule.fitness`. The :class:`Normalization` instance
does this automatically.

Normalization functions calculate the fitness value of a molecule and
place it in :attr:`.MacroMolecule.fitness`. Multiple normalization
functions can be applied in sequence. Only the last normalization
function applied needs to place a value between 0 (exlusive) and
infinity into :attr:`.MacroMolecule.fitness`. The others can do
whatever scaling and data manipulation is necessary.

If a normalization function does not fit neatly into a single function
make sure that any helper functions are private, i.e. that their names
start with a leading underscore.

"""

import numpy as np
import copy
import logging

from .ga_population import GAPopulation


logger = logging.getLogger(__name__)


class Normalization:
    """
    A class for carrying out normalization of fitness values.

    Attributes
    ----------
    funcs : :class:`list` of :class:`.FunctionData`
        Holds all the normalization functions to be applied each
        generation, in the order in which they are to be applied.

    """

    def __init__(self, funcs):
        """
        Initializes a :class:`Normalization` instance.

        Parameters
        ----------
        funcs : :class:`list` of :class:`.FunctionData`
            Holds all the normalization functions to be applied each
            generation, in the order in which they are to be applied.

        """

        self.funcs = funcs

    def __call__(self, population):
        """
        Applies the normalization functions on `population`.

        Parameters
        ----------
        population : :class:`.GAPopulation`
            The population whose members need to have their fitness
            values normalized.

        Returns
        -------
        None : :class:`NoneType`

        """

        fitness_func = population.ga_tools.fitness.name
        # First make sure that all the fitness values are reset and
        # hold the value of the approraite fitness function.
        for macro_mol in population:
            # If the unscaled fitness value is None it means that the
            # fitness calculation failed - assign the minimum fitness
            # value.
            if macro_mol.unscaled_fitness[fitness_func] is None:
                macro_mol.fitness = 1e-4
            else:
                macro_mol.fitness = copy.deepcopy(
                            macro_mol.unscaled_fitness[fitness_func])

        # Make a population of members where all fitness values are
        # valid. No point in normalizing molecules whose fitness value
        # was ``None``. The advantage of this is that when writing
        # normalization functions you can assume all fitness values
        # have the same type.

        valid_pop = GAPopulation(*(mol for mol in population if
                                 mol.unscaled_fitness[fitness_func]
                                 is not None))

        # Remove duplicates, otherwise the normalization function may
        # be applied to the same molecule twice in a single step.
        valid_pop.remove_duplicates()

        # If there were no valid molecules, no need to normalize.
        if len(valid_pop) == 0:
            return

        for func_data in self.funcs:
            logger.debug('Applying "{}()".'.format(func_data.name))
            getattr(self, func_data.name)(valid_pop,
                                          **func_data.params)

    def cage(self, population, cavity, window):
        """
        A normalization function for use with :func:`.cage`

        Parameters
        ----------
        population : :class:`.GAPopulation`
            The population whose members need to have their fitness
            values normalized.

        cavity : :class:`float`
            The desired size of the cage cavity.

        window : :class:`float`
            The desired size of the largest cage window.

        Returns
        -------
        None : :class:`NoneType`

        """

        for mem in population:
            cavity_diff = abs(mem.fitness[0] - cavity)
            window_diff = abs(mem.fitness[1] - window)
            mem.fitness = [cavity_diff, window_diff,
                           *mem.fitness[2:]]

    def combine(self, population, coefficients, exponents):
        """
        Combines elements in the :attr:`.MacroMolecule.fitness`.

        This function assumes that :attr:`.MacroMolecule.fitness`
        of the population's members is an array. It raises the elements
        of the array to the powers in `exponents` and then multiplies
        them by the values in `coefficients`. Lastly, the individual
        array elements are summed to create a final value.

        Parameters
        ----------
        population : :class:`.GAPopulation`
            The population whose members need to have their fitness
            values normalized.

        coefficients : :class:`list` of :class:`int` or :class:`float`
            Before summing all the elements in
            :attr:`~.MacroMolecule.fitness`, they are multiplied by the
            values in this list.

        exponents : :class:`list` of :class:`int` or :class:`float`
            Before summing all the elements in the
            :attr:`~.MacroMolecule.fitness`, they are raised to the
            powers in this list.

        Returns
        -------
        None : :class:`NoneType`

        """

        for macro_mol in population:
            new_array = np.power(macro_mol.fitness, exponents)
            new_array = np.multiply(new_array, coefficients)
            macro_mol.fitness = sum(new_array)

    def divide(self, population, val):
        """
        Divides each fitness value by `val`.

        Parameters
        ----------
        population : :class:`.GAPopulation`
            The population being normalized.

        val : :class:`int` or :class:`float` or :class:`numpy.ndarray`
            The value by which each fitness value is divided.

        Returns
        -------
        None : :class:`NoneType`

        """

        for mem in population:
            mem.fitness /= val

    def magnitudes(self, population):
        """
        Sacles values of elements in :attr:`.MacroMolecule.fitness`.

        This normalization function assumes that the value in
        :attr:`~.MacroMolecule.fitness` is a numpy array:

        .. code-block:: python

            macro_mol.fitness = np.array([1000, 5e-5, 25])

        Notice that each element has a very different order of
        magnitude. For example, if the fitness function calculated the
        energy value of the molecule, its radius and the number of
        atoms in it, this could be the data placed into
        :attr:`~.MacroMolecule.fitness` as an array.

        When calculating the total fitness based on these values, it
        would be useful to make them comparable. As it is, you can't
        compare 10,000 kJ mol-1 and 5e-5 Angstroms. However, what you
        can do is check how much bigger or smaller than the average
        10,000 kJ mol-1 and 5e-5 Angstrom are. Then replace these
        values with the size divided by the average. For example:

        .. code-block:: python

            average_vals = np.array([500, 1e-5, 5])
            macro_mol.fitness = macro_mol.fitness / average_vals
            macro_mol.fitness  # np.array([2, 0.5, 5])

        It shows that the energy of a given `macro_mol` is twice as
        large as the mean energy of the population. The radius is half
        the average molecular radius of the population and so on.

        Now these values can be combined in a reasonable way. However,
        that will have to be done by other normalization functions.
        This one only scales relative to the population average.

        Parameters
        ----------
        population : :class:`.GAPopulation`
            The population whose fitness values are normalized.

        Returns
        -------
        None : :class:`NoneType`

        """

        # Get the mean of each element.
        means = population.mean(lambda x: x.fitness)
        logger.debug('Means used in magnitudes: {}'.format(means))

        for macro_mol in population:
            macro_mol.fitness = macro_mol.fitness / means

    def shift_elements(self, population, indices):
        """
        Makes all elements in :attr:`.MacroMolecule.fitness` positive.

        Assume you have a fitness array,

        .. code-block:: python

            macro_mol.fitness = np.array([1, -10, 1])

        One way to convert the fitness array into a fitness value is
        by summing the elements (see :meth:`Normalization.combine` for
        this). The result would be:

        .. code-block:: python

            # After summing all values.
            macro_mol.fitness  # -8

        Clearly this doesn't work, because the resulting fitness value
        is not a positive number. To fix this, the ``-10`` should be
        shifted to a positive value.

        This normalization function looks at the elements specified by
        by `indices`. It then finds the minimum value of these elements
        in the population. It then shifts the elements by this value.

        For example, take a population with the fitness values

        .. code-block:: python

            mol1.fitness = np.array([1, -5, 5])
            mol2.fitness = np.array([3, -10, 2])
            mol3.fitness = np.array([2, 20, 1])

        If the value of `indices` was ``[1]`` then after this
        normalization function was applied the result would be

        .. code-block:: python

            mol1.fitness  # np.array([1, 5.1, 5])
            mol2.fitness  # np.array([3, 0.1, 2])
            mol3.fitness  # np.array([2, 30.1, 1])

        Notice that all in all the fitness arrays the second element
        has been shifted by ``-10.1``. This is the smallest value of the
        second element in the population multiplied by ``1.01`` to
        prevent any of the values being ``0``.

        If `indices` was ``[0,1]``

        .. code-block:: python

            mol1.fitness  # np.array([2.01, 5.1, 5])
            mol2.fitness  # np.array([4.01, 0.1, 2])
            mol3.fitness  # np.array([3.01, 30.1, 1])

        The result is the same as before, only now all the first and
        second elements have been shifted.

        Parameters
        ----------
        indices : :class:`list` of :class:`int`
            This holds the indices of elements in the
            :attr:`~.MacroMolecule.fitness` array which should be
            shifted.

        Returns
        -------
        None : :class:`NoneType`

        """

        # Get all the fitness arrays a matrix.
        fmat = np.array([x.fitness for x in population])

        # Get the minimum values of each element in the population.
        mins = np.min(fmat, axis=0)
        # Convert all the ones which are not to be shifted to 0 and
        # multiply the which are to be shifted by 1.01.
        shift = np.zeros(len(mins))
        for ind in indices:
            shift[ind] = 1.01
        shift = abs(np.multiply(mins, shift))

        for macro_mol in population:
            macro_mol.fitness += shift

    def invert(self, population):
        """
        Raises fitness values to the power of ``-1``.

        Parameters
        ----------
        population : :class:`.GAPopulation`
            The population to be normalized.

        Returns
        -------
        None : :class:`NoneType`

        """

        for macro_mol in population:
            macro_mol.fitness = 1 / macro_mol.fitness
