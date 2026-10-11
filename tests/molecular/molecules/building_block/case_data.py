from collections.abc import Sequence
from dataclasses import dataclass

import stk


@dataclass(slots=True, frozen=True)
class CaseData:
    building_block: stk.BuildingBlock
    functional_groups: Sequence[stk.FunctionalGroup]
    known_repr: str
    core_atom_ids: Sequence[int]
    placer_ids: Sequence[int]
