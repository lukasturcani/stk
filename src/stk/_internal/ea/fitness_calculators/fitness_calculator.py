import typing

from stk._internal.ea.molecule_records.molecule import MoleculeRecord


class FitnessCalculator:
    """
    Abstract base class for fitness value calculators.

    Examples:

        *Subclass Implementation*

        You only need to implement :meth:`.get_fitness_value`. The source
        cod of any of the classes listed in :mod:`.fitness_calculator`, can
        serve as good examples.
    """

    def get_fitness_value(self, record: MoleculeRecord) -> typing.Any:
        """
        Return the fitness value of a molecule.

        Parameters:
            record:
                The molecule whose fitness value should be calculated.

        Returns:
            The fitness value.
        """
        raise NotImplementedError()
