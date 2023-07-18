import pathlib
import typing

import stk

from .case_data import CaseData

T = typing.TypeVar("T", bound=stk.MoleculeRecord)


def test_selection_plotter(
    tmp_path: pathlib.Path,
    case_data: CaseData,
) -> None:
    _test_selection_plotter(
        selector=case_data.selector,
        population=case_data.population,
        filename=tmp_path / "selection",
    )


def _test_selection_plotter(
    selector: stk.Selector[T],
    population: dict[T, float],
    filename: pathlib.Path,
) -> None:
    stk.SelectionPlotter(filename, selector)
    tuple(selector.select(population))
