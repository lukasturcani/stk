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
import sys

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

        for func_data in self.funcs:
            scaling_func = partial(getattr(self, func_data.name),
                                    **func_data.params)
            scaling_func(population)

    def combine(self, population, elements, coefficients, exponents):
        """
        Combines elements in the unscaled_fitness attribute of members.

        This normalization function is not stand-alone. It is desinged
        to be chained to other normalization functions.

        This function assumes that the `unscaled_fitness` attribute
        of the population's members is an array. It multiplies the
        members of the array by the values in `coefficients` and then
        raises them to the values in `exponents`. Lastly, the
        individual array elements are summed to create a new smaller
        array.

        Parameters
        ----------
        population : Population
            The population whose members need to have their fitness
            values normalized.

        elements : list of ints
            If the list is [3, 2, 3] then the first 3 elements are
            summed to for the first element of the normalized vector.
            The second element of the normalized vector is formed by
            summing the next 2 elements of the original. Finally, the
            last element of the normalized array is made by summing the
            next 3 elements.

        coefficients : list of ints or floats
            Before summing all the elements in the `unscaled_fitness`
            attribute of members, their are multiplied by the numbers
            in this list.

        exponents : list of ints or floats
            Before summing all the elements in the `unscaled_fitness`
            attribute of members, their are raised to the numbers in
            this list.

        Modifies
        --------
        fitness : numpy.array
            This attribute is altered for the populations members.

        Returns
        -------
        None : NoneType

        """

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

        When calculating the total fitness pased on thse values, it
        would be useful to make them comparable. As it is, you can't
        compare 10,000 kJ mol-1 and 5e-5 Angstroms. However what you
        can do is check how much bigger or smaller than average
        10,000 kJ mol-1 and 5e-5 Angstrom are. Then replace these
        values with the size relative to the average. For example:

            macro_mol.fitness = np.array([2, 0.5, 3])

        This would be the outout of this function. It shows that the
        energy of a given `macro_mol` is twice as large as the mean
        energy of the population. The radius is half the average
        molecular radius of the population and so on.

        Now these values can be combined in a reasonable way. However,
        that will have to be done by other normalization functions.
        This one only scales relative to the population average.

        Assuming that the populations members have a numpy array in
        their `fitness` attribute. The following steps are performed:

            1) Calculate the mean value for each element across the
               population.

            2) Replace each element with its difference from the mean.

            3) Calculate the standard deviation of each element across
               the population.

            4) Replace each element with its value divided by the
               standard deviation.


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
        means = population.mean(lambda x : x.fitness)

        # Replace values by deviations from mean.
        for macro_mol in population:
            macro_mol.fitness -= means

        # Calculate the standard devation of each element.

        # First, get a matrix where each row consistents of the fitness
        # parameters of a population member.
        pop_mat = np.array([x.unscaled_fitness for x in population])

        # To get the standard deviation, square each element and then
        # sum all the squares of elements in the same column. Divide
        # by the number of rows and square root.
        devs = np.sqrt(np.sum(np.square(pop_mat), axis=0) / len(pop))

        # Update the fitness values.
        for macro_mol in population:
            macro_mol.fitness = macro_mol.fitness / devs

    def pareto(self, population):
        """

        """

        ...
