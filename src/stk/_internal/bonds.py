from dataclasses import dataclass

from stk._internal.atom import Atom


@dataclass(frozen=True, slots=True)
class IntegerBond:
    atom1: Atom
    atom2: Atom
    order: int


@dataclass(frozen=True, slots=True)
class PeriodicBond:
    atom1: Atom
    atom2: Atom
    order: int
    periodicity: tuple[int, int, int]


@dataclass(frozen=True, slots=True)
class DativeBond:
    atom1: Atom
    atom2: Atom
    order: int
