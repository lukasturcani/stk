"""
Fitness Calculator
==================

"""


class FitnessCalculator:
    """
    Calculates fitness values of molecules.

    """

    def get_fitness_value(self, molecule):
        """
        Return the fitness value of `molecule`.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule whose fitness value should be calculated.

        Returns
        -------
        :class:`object`
            The fitness value of `molecule`.

        """

        raise NotImplementedError()
