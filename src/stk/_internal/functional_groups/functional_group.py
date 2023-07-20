from collections.abc import Sequence
from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class FunctionalGroup:
    bonders: Sequence[int]
    deleters: Sequence[int]
