import itertools
import typing
from collections.abc import Sequence

import stk

from .case_data import CaseData

T = typing.TypeVar("T", bound=stk.MoleculeRecord)


def test_select(case_data: CaseData) -> None:
    _test_select(
        selector=case_data.selector,
        population=case_data.population,
        selected=case_data.selected,
    )


def _test_select(
    selector: stk.Selector[T],
    population: dict[T, float],
    selected: Sequence[stk.Batch[T]],
) -> None:
    inchi = stk.Inchi()
    for batch1, batch2 in itertools.zip_longest(
        sorted(selector.select(population), key=get_inchis),
        sorted(selected, key=get_inchis),
    ):
        inchis1 = tuple(
            inchi.get_key(record.get_molecule()) for record in batch1
        )
        inchis2 = tuple(
            inchi.get_key(record.get_molecule()) for record in batch2
        )
        assert inchis1 == inchis2


def get_inchis(batch: stk.Batch[T]) -> tuple[str, ...]:
    inchi = stk.Inchi()
    return tuple(inchi.get_key(record.get_molecule()) for record in batch)
