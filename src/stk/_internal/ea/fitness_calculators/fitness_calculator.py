import typing

from stk._internal.ea.molecule_record import MoleculeRecord

T = typing.TypeVar("T", bound=MoleculeRecord)


class FitnessCalculator(typing.Generic[T]):
    """
    Abstract base class for fitness value calculators.

    Examples:

        *Subclass Implementation*

        You only need to implement :meth:`.get_fitness_value`.
    """

    def get_fitness_value(self, record: T) -> typing.Any:
        """
        Return the fitness value of a molecule.

        Parameters:
            record:
                The molecule whose fitness value should be calculated.
        Returns:
            The fitness value.
        """
        raise NotImplementedError()
