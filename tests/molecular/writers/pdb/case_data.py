import stk

from dataclasses import dataclass


@dataclass(slots=True, frozen=True)
class CaseData:
    """
    A :class:`.PdbWriter` test case.

    """

    molecule: stk.Molecule
    writer: stk.PdbWriter
    string: str
    periodic_info: stk.PeriodicInfo | None
