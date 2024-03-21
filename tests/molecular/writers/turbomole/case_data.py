import stk

from dataclasses import dataclass


@dataclass(slots=True, frozen=True)
class CaseData:
    """
    A :class:`.TurbomoleWriter` test case.

    """

    molecule: stk.Molecule
    writer: stk.TurbomoleWriter
    string: str
    periodic_info: stk.PeriodicInfo | None
