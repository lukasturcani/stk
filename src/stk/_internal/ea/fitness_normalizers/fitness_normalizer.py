import typing
from typing import Any

T = typing.TypeVar("T")


class FitnessNormalizer(typing.Generic[T]):
    """
    Abstract base class for fitness normalizers.

    A fitness normalizer takes a dictionary mapping
    :class:`.MoleculeRecord` instances to their fitness values
    and returns a new dictionary mapping them to their normalized
    fitness values. The primary benefit of a normalizer vs a
    :class:`.FitnessCalculator` is that a :class:`.FitnessNormalizer`
    has access to all members in the population when it is calculating
    the normalized fitness value, whereas a :class:`.FitnessCalculator`
    does not.
    """

    def normalize(self, fitness_values: dict[T, Any]) -> dict[T, Any]:
        """
        Normalize some fitness values.

        Parameters:
            fitness_values:
                The molecules which need to have their fitness values
                normalized.
        Returns:
            The new fitness value for each molecule.
        """
        raise NotImplementedError()
