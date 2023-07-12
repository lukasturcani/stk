import pathlib

import pandas as pd
import stk

from .case_data import CaseData


def test_progress_plotter(tmp_path: pathlib.Path, case_data: CaseData) -> None:
    _test_progress_plotter(
        plotter=case_data.plotter,
        path=tmp_path / "plot.png",
        plot_data=case_data.plot_data,
    )


def _test_progress_plotter(
    plotter: stk.ProgressPlotter,
    path: pathlib.Path,
    plot_data: pd.DataFrame,
) -> None:
    plotter.write(path)
    assert plotter.get_plot_data().equals(plot_data)
