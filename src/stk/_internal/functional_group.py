from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class FunctionalGroup:
    bonders: list[int]
    deleters: list[int]
