import typing
from dataclasses import dataclass

import stk

T = typing.TypeVar("T", bound=stk.MoleculeRecord)


@dataclass(frozen=True, slots=True)
class CaseData(typing.Generic[T]):
    selector: stk.Selector[T]
    population: dict[T, float]
