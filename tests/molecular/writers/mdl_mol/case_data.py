import stk

from dataclasses import dataclass


@dataclass(slots=True, frozen=True)
class CaseData:
    molecule: stk.Molecule
    writer: stk.MolWriter
    string: str
