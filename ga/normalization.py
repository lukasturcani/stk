"""
Defines normalization functions via the Normalization class.

Normalization functions are functions that recalculate the fitness
values of members in a population. The difference between fitness
and normalization functions is that fitness functions are only use the
MacroMolecule to calculate its fitness. Normalzation functions have
access to all MacroMolecules in a population however. As a result they
can scale the fitness values across the entire population. This is
useful if you want to ensure a spread of fitness values in the
population.

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

Normalization functions calculate the fitness value of a molecule and
place it in the `fitness` attribute. This is not universally true, if a
function is intended to be ``chained`` to other normalization
functions. In other words, if the normalization function is meant to be
applied before another one is. This is done via the ``chain()`` method.

Note that fitness functions are only executed once per MacroMolecule,
while normalization functions can be executed multiple times. Because
fitness functions are only executed once, the value they calculate
should never be overwritten. As a result when using normalization
functions, the fitness functions should write to the `unscaled_fitness`
attribute of MacroMolecules. This will be unchanged once calculated.
The normalization functions can then use the value in
`unscaled_fitness` to calculate a value for `fitness`. The calculated
`fitness` value may change each generation, depending on the scaling
procedure used.

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

    def __init__(self, func_data):
        """
        Initializes a Normalization instance.

        Parameters
        ----------
        func_data : FunctionData
            A FunctionData

        """

        self.scaling_func = partial(getattr(self, func_data.name),
                                    **func_data.params)

    def __call__(self, population):
        """
        Applies the normalization function on `population`.

        Parameters
        ----------
        population : Population
            The population whose members need to have their fitness
            values normalized.

        """

        self.scaling_func(population)

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
        unscaled_fitness : numpy.array
            This attribute is altered for the populations members.

        Returns
        -------
        None : NoneType

        """



    @staticmethod
    def carrots_and_sticks(population, carrot_coeffs, stick_coeffs,
                           carrot_exponents, stick_exponents):
        """
        Applies the ``carrots and sticks`` normalization.

        This function requires that the fitness functions places a tuple
        of 2 arrays in the `unscaled_fitness` attribute of members.

            mem.unscaled_fitness =  (carrots, sticks)

        where

            carrots = np.array([c1, c2, c3])

        and

            sticks = np.array([s1, s2])

        Note that the arrays can be of any size.

        The goal is to maximize the carrot values and minimize the
        stick values. As a result, the fitness of an individual is given
        by

            (1) fitness = sum(carrots) + 1/sum(sticks)

        What if c1 is 1000 and c2 is 0.01?

        This means that the fitness value is dominated by c1 and c2
        is barely getting optimized. To fix this the values of all
        elements are recalculated

            (2) Se = e / <e>

        where ``e`` can represent a carrot or stick parameter (c1, s2
        etc.). The <e> is the average of that parameter across all
        members in the population.

        This scaling means that c1 and c2 are rescaled to be around the
        same order of magnitude.

        What if we want c1 to twice as important to fitness as c3?

        After equation (2) the elements are rescaled again

            (3) Ne = A*(Se^a)

        Here the ``Ne`` represents the parameter ``e`` after the second
        rescaling.

        This means that the sums showin in equation (1) are

            (4) sum(carrots) = A*(Nc1^a) + B*(Nc2^b) + C*(Nc3^b)

        The prefix ``N`` in front of c1, c1 and c3 means that the
        operations in equations (2) and (3) have been applied to the
        parameters.

        So if you want c1 to be twice as important to fitness as c2
        set

            A = 1 and B = 2

        This is done via the `coeff` and `exponents` parameters.

        Parameters
        ----------
        population : Population
            The population whose fitness values are normalized.

        carrot_coeffs : numpy.array
            The coeffients of the carrot parameters.

        stick_coeffs : numpy.array
            The coefficients of the stick parameters.

        carrot_exponents : numpy.array
            The exponents of the carrot parameters.

        stick_exponents : numpy.array
            The exponents of the stick parameters.

        Modifies
        --------
        fitness : float
            This attribute in all of the population's members is
            modified.

        Returns
        -------
        None : NoneType

        """

        unscaled_carrots = [x.unscaled_fitness[0] for x in population if
                            isinstance(x.unscaled_fitness, tuple)]

        unscaled_sticks = [x.unscaled_fitness[1] for x in population if
                           isinstance(x.unscaled_fitness, tuple)]

        _carrot_means = np.mean(unscaled_carrots, axis=0)
        if not isinstance(_carrot_means, np.ndarray):
            raise TypeError(('Unscaled fitness values do not have'
                             ' appropriate type.'))

        carrot_means = []
        for x in _carrot_means:
            if x == 0:
                carrot_means.append(1)
            else:
                carrot_means.append(x)

        _stick_means = np.mean(unscaled_sticks, axis=0)
        if not isinstance(_carrot_means, np.ndarray):
            raise TypeError(('Unscaled fitness values do not have'
                             ' appropriate type.'))

        stick_means = []
        for x in _stick_means:
            if x == 0:
                stick_means.append(1)
            else:
                stick_means.append(x)

        for macro_mol in population:

            # If one or more of the fitness parameters failed,
            # return minimum fitness.
            if macro_mol.fitness_fail:
                macro_mol.fitness = 1e-4
                continue

            # Calculate the scaled fitness parameters by dividing the
            # unscaled ones by the fitness.
            scaled_carrots = np.divide(macro_mol.unscaled_fitness[0],
                                       carrot_means)

            scaled_sticks = np.divide(macro_mol.unscaled_fitness[1],
                                      stick_means)

            try:
                scaled_carrots = np.power(scaled_carrots,
                                          carrot_exponents)
                scaled_sticks = np.power(scaled_sticks,
                                         stick_exponents)

                scaled_carrots = np.multiply(scaled_carrots,
                                             carrot_coeffs)
                scaled_sticks = np.multiply(scaled_sticks,
                                            stick_coeffs)

            # If the user forgot put the wrong number values in the
            # exponent or coefficient arrays MMEA will tell them and
            # exit gracefully.
            except  ValueError:
                print(('Fitness function calculates carrot array of '
                       'size {} and stick array of size {}. This does '
                       'not match the array sizes given to the '
                       'normalization function:\n'
                       '\tcarrot_coeffs : {}\n'
                       '\tcarrot_exponents : {}\n'
                       '\tstick_coeffs : {}\n'
                       '\tstick_exponents : {}\n').format(
                                   len(scaled_carrots),
                                   len(scaled_sticks),
                                   len(carrot_coeffs),
                                   len(carrot_exponents),
                                   len(stick_coeffs),
                                   len(stick_exponents)))
                sys.exit()

            carrot_term = np.sum(scaled_carrots)
            penalty_term = np.sum(scaled_sticks)
            penalty_term =  np.divide(1,penalty_term)
            if penalty_term > 1e101:
                penalty_term = 1e101

            macro_mol.fitness = penalty_term + carrot_term


    def chain(self, population, funcs):
        """
        Applies multiple normalization functions, one after the other.

        Parameters
        ----------
        population : Population
            The population whose fitness values are normalized.

        funcs : list of FunctionData objects
            The FunctionData objects representing normalization
            functions to be applied, in the order in which they should
            be applied.

        Returns
        -------
        None : NoneType

        """

        for func_data in funcs:
            func = getattr(self, func_data.name)
            func(population, **func_data.params)

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
        unscaled_fitness : numpy.array
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
