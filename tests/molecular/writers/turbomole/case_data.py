import stk

from dataclasses import dataclass


@dataclass(slots=True, frozen=True)
class CaseData:
    molecule: stk.Molecule
    writer: stk.TurbomoleWriter
    string: str
    periodic_info: stk.PeriodicInfo | None
