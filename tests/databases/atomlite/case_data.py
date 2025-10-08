from dataclasses import dataclass

import stk


@dataclass(slots=True, frozen=True)
class CaseData:
    molecules: list[stk.Molecule]
    property_dicts: list[dict]
    expected_count: int
    name: str
