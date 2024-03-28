import stk
from dataclasses import dataclass
from collections.abc import Sequence


@dataclass(slots=True, frozen=True)
class CaseData:
    building_block: stk.BuildingBlock
    functional_groups: Sequence[stk.FunctionalGroup]
    known_repr: str
    core_atom_ids: Sequence[int]
    placer_ids: Sequence[int]
