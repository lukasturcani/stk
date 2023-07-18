import typing
from collections.abc import Sequence
from dataclasses import dataclass

import stk

T = typing.TypeVar("T", bound=stk.MoleculeRecord)


@dataclass(frozen=True, slots=True)
class CaseData(typing.Generic[T]):
    selector: stk.Selector[T]
    population: dict[T, float]
    selected: Sequence[stk.Batch[T]]

    @staticmethod
    def new(
        selector: stk.Selector[T],
        population: dict[T, float],
        selected: Sequence[Sequence[int]],
    ) -> "CaseData":
        records = list(population)
        batches = []
        for batch in selected:
            batches.append(
                stk.Batch(
                    records={
                        records[member]: population[records[member]]
                        for member in batch
                    },
                    key_maker=stk.Inchi(),
                )
            )
        return CaseData(selector, population, batches)
