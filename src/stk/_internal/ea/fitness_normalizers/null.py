import typing
from typing import Any

from .fitness_normalizer import FitnessNormalizer

T = typing.TypeVar("T")


class NullFitnessNormalizer(FitnessNormalizer[T]):
    """
    Does nothing.

    This normalizer just yields the molecule records passed to it,
    without changing them in any way.
    """

    def normalize(self, fitness_values: dict[T, Any]) -> dict[T, Any]:
        return fitness_values
