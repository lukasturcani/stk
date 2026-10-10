from dataclasses import dataclass

import stk


@dataclass(slots=True, frozen=True)
class CaseData:
    molecule: stk.Molecule
    writer: stk.PdbWriter
    string: str
    periodic_info: stk.PeriodicInfo | None
