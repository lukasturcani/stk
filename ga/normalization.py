"""
Defines normalization functions via the Normalization class.

Normalization functions are functions that recalculate the fitness
values of members in a population. The difference between fitness
and normalization functions is that fitness functions are only use the
MacroMolecule to calculate its fitness. Normalzation functions have
access to all MacroMolecules in a population. As a result they can
scale the fitness values across the entire population. This is useful
if you want to ensure a spread of fitness values in the population.

Each generation, a number of normalizatoin functions can be applied
in sequence.

Extending MMEA: Adding normalization functions.
-----------------------------------------------
If a new normalization function is to be added to MMEA it should be
added as a method in the ``Normalization`` class defined in this
module. The only requirements are that the first argument is
``population`` (excluding ``self``).

The naming requirement exists to help users identify which arguments
are handled automatically by MMEA and which they need to define in the
input file. The convention is that if the normalization function takes
an argument called  ``population`` it does not have to be specified in
the input file.

Normalization functions should not interact with the `unscaled_fitness`
attribute in any way. They should only modify the value in `fitness`.
Before the first normalization function is applied each generation, the
value in `unscaled_fitness` is copied into `fitness`. This happens
automatically.

Normalization functions calculate the fitness value of a molecule and
place it in the `fitness` attribute. Multiple normalization functions
can be applied in sequence. Only the last normalization function
applied needs to place a value between 0 (exlusive) and infinity in the
`fitness` attribute. The others can do whatever scaling is necessary
to for the problem at hand.

If a normalization function does not fit neatly into a single function
make sure that any helper functions are private, ie that their names
start with a leading underscore.

"""


from functools import partial
import numpy as np
import sys, copy, logging

from .population import Population


logger = logging.getLogger(__name__)


class Normalization:
    """
    A class for carrying out normalization of fitness values.

    Attributes
    ----------
    scaling_func : functools. partial

    """

    def __init__(self, funcs):
        """
        Initializes a Normalization instance.

        Parameters
        ----------
        funcs : list of FunctionData instances
            Holds all the normalization functions to be applied each
            generation, in the order in which they are to be applied.

        """

        self.funcs = funcs

    def __call__(self, population):
        """
        Applies the normalization function on `population`.

        Parameters
        ----------
        population : Population
            The population whose members need to have their fitness
            values normalized.

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

        valid_pop = Population(*(mol for mol in population if
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
        A normalization function for the `cage` fitness function.

        Parameters
        ----------
        population : Population
            The population whose members need to have their fitness
            values normalized.

        cavity : float
            The desired size of the cage cavity.

        window : float
            The desired size of the largest cage window.

        Modifies
        --------
        fitness : numpy.array
            The fitness attribute of members is updated.

        Returns
        -------
        None : NoneType

        """

        for mem in population:
            cavity_diff = abs(mem.fitness[0] - cavity)
            window_diff = abs(mem.fitness[1] - window)
            mem.fitness = [cavity_diff, window_diff,
                           *mem.fitness[2:]]

    def combine(self, population, coefficients, exponents):
        """
        Combines elements in the `fitness` attribute of members.

        This function assumes that the `fitness` attribute
        of the population's members is an array. It raises the members
        of the array to the values in `exponents` and then multiplies
        them with the values in `coefficients`. Lastly, the individual
        array elements are summed to create a final value.

        Parameters
        ----------
        population : Population
            The population whose members need to have their fitness
            values normalized.

        coefficients : list of ints or floats
            Before summing all the elements in the `fitness` attribute,
            their values are multiplied by the numbers in this list.

        exponents : list of ints or floats
            Before summing all the elements in the `fitness` attribute,
            their values are raised to the numbers in this list.

        Modifies
        --------
        fitness : numpy.array
            This attribute is altered for the populations members. It
            is changed from an array to a float.

        Returns
        -------
        None : NoneType

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
        population : Population
            The population being normalized.

        val : numerical or numpy.array of numericals
            The value by which each fitness value is divided.

        Modifies
        --------
        fitness : numerical or numpy.array of numericals
            Each fitness value of a population member is divided by
            `val`.

        Returns
        -------
        None : NoneType

        """

        for mem in population:
            mem.fitness /= val

    def magnitudes(self, population):
        """
        Normalizes the relative values of elements in `fitness`.

        This normalization function assumes that the value in `fitness`
        is a numpy array:

            macro_mol.fitness = np.array([1000, 5e-5, 25])

        Notice that each element has a very different order of
        magnitude. For example if the fitness function calculated the
        energy value of the molecule, its radius and the number of
        atoms in it, this could be the data placed in the `fitness`
        attribute.

        When calculating the total fitness based on these values, it
        would be useful to make them comparable. As it is, you can't
        compare 10,000 kJ mol-1 and 5e-5 Angstroms. However what you
        can do is check how much bigger or smaller than average
        10,000 kJ mol-1 and 5e-5 Angstrom are. Then replace these
        values with the size relative to the average. For example:

            macro_mol.fitness = np.array([2, 0.5, 3])

        This would be the output of this function. It shows that the
        energy of a given `macro_mol` is twice as large as the mean
        energy of the population. The radius is half the average
        molecular radius of the population and so on.

        Now these values can be combined in a reasonable way. However,
        that will have to be done by other normalization functions.
        This one only scales relative to the population average.

        Parameters
        ----------
        population : Population
            The population whose fitness values are normalized.

        Modifies
        --------
        fitness : numpy.array
            This attribute is altered for the populations members.

        Returns
        -------
        None : NoneType

        """

        # Get the mean of each element.
        means = population.mean(lambda x: x.fitness)
        logger.debug('Means used in magnitudes: {}'.format(means))

        for macro_mol in population:
            macro_mol.fitness = macro_mol.fitness / means

    def shift_elements(self, population, indices):
        """
        Maps elements in `fitness` array to positive values.

        Assumy you have a fitness array,

            macro_mol.fitness = [1, -10, 1]

        One way to convert the fitness array into a fitness value is
        by summing the elements (see the combine() normalization
        function for this)

            macro_mol.fitness = -8

        Clearly this doesn't work, because the resulting fitness value
        is not a positive number. To fix this the -10 should be shifted
        to a positive value.

        This normalization function looks at the elements specified by
        by `indices`. It then finds the minimum value of these elements
        in the population. It then shifts the elements by this value.

        For example, take a population of

            mol1.fitness = [1, -5, 5]
            mol2.fitness = [3, -10, 2]
            mol3.fitness = [2, 20, 1]

        If the value of `indices` was [1] the after this normalization
        function was applied the result would be

            mol1.fitness = [1, 5.1, 5]
            mol2.fitness = [3, 0.1, 2]
            mol3.fitness = [2, 30.1, 1]

        If `indices` was [0,1]

            mol1.fitness = [2.01, 5.1, 5]
            mol2.fitness = [4.01, 0.1, 2]
            mol3.fitness = [3.01, 30.1, 1]

        The shift is  1.01 * magnitude of smallest value. This prevents
        0s in the array.

        Parameters
        ----------
        indices : list of ints
            This holds the indices of elements in the `fitness`
            array which should be shifted to positive values.

        Modifies
        --------
        fitness : numpy.array
            The `fitness` atttribute of the population's members is
            changed in accordance to the docstring.

        Returns
        -------
        None : NoneType

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
        Convertes a fitness value to 1/fitness.

        Parameters
        ----------
        population : Population
            The population to be normalized.

        Modifies
        --------
        fitness : float
            The `fitness` attribute of the population's members is
            changed to 1/`fitness`.

        Returns
        -------
        None : NoneType

        """

        for macro_mol in population:
            macro_mol.fitness = 1 / macro_mol.fitness
