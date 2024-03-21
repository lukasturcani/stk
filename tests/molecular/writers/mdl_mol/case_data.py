import stk

from dataclasses import dataclass


@dataclass(slots=True, frozen=True)
class CaseData:
    """
    A :class:`.MolWriter` test case.

    """

    molecule: stk.Molecule
    writer: stk.MolWriter
    string: str
