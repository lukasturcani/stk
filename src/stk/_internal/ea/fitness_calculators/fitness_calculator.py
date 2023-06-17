"""
Fitness Calculator
==================

.. toctree::
    :maxdepth: 2

    Fitness Function <stk.ea.fitness_calculators.fitness_function>
    Property Vector <stk.ea.fitness_calculators.property_vector>

"""


class FitnessCalculator:
    """
    Abstract base class for fitness value calculators.

    Examples
    --------
    *Subclass Implementation*

    You only need to implement :meth:`.get_fitness_value`. The source
    cod of any of the classes listed in :mod:`.fitness_calculator`, can
    serve as good examples.

    """

    def get_fitness_value(self, molecule):
        """
        Return the fitness value of `molecule`.

        Parameters
        ----------
        molecule : :class:`.ConstructedMolecule`
            The molecule whose fitness value should be calculated.

        Returns
        -------
        :class:`object`
            The fitness value of `molecule`.

        """

        raise NotImplementedError()
