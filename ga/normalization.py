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
        raises them to the values in `exponents`. Lastly, the individual
        array elements are summed to create a new smaller array.

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

    def magnitudes(population, use_fitness=False):
        """
        Normalizes the values in unscaled_fitness by the average.

        This normalization function is not stand-alone. It is desinged
        to be chained to other normalization functions.

        Assumes that the populations members have a numpy array in their
        `unscaled_fitness` attribute. The following steps are performed:

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

        use_fitness : bool (default = False)
            If ``True`` the normalization function is applied to the
            values in the `fitness` attribute of members instead of
            `unscaled_fitness`. This is necessary when using multiple
            normalization functions in a row.

        Modifies
        --------
        fitness : numpy.array
            This attribute is altered for the populations members.

        Returns
        -------
        None : NoneType

        """

        # Ge the name of the attribute which is
        if use_fitness:
            attr_name

        # Get the mean of each element.
        means = population.mean(lambda x : x.unscaled_fitness)
