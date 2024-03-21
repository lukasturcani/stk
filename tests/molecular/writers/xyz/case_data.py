import stk

from dataclasses import dataclass


@dataclass(slots=True, frozen=True)
class CaseData:
    """
    A :class:`.XyzWriter` test case.

    """

    molecule: stk.Molecule
    writer: stk.XyzWriter
    string: str
