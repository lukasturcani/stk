import stk
from dataclasses import dataclass


@dataclass(slots=True, frozen=True)
class CaseData:
    """
    A test case.

    """

    building_block: stk.BuildingBlock
    functional_groups: tuple
    known_repr: str
    core_atom_ids: tuple[int, ...]
    placer_ids: tuple[int, ...]
