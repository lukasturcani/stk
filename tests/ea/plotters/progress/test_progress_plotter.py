def test_progress_plotter(tmp_path, case_data):
    """
    Test :class:`.ProgressPlotter`.

    Parameters
    ----------
    tmp_path : :class:`pathlib2.Path`
        A directory into which the plot is written.

    case_data : :class:`.CaseData`
        A test case. Holds the plotter to test and the correct plot
        data.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_progress_plotter(
        plotter=case_data.plotter,
        path=tmp_path / "plot.png",
        plot_data=case_data.plot_data,
    )


def _test_progress_plotter(plotter, path, plot_data):
    """
    Test :class:`.ProgressPlotter`.

    Parameters
    ----------
    plotter : :class:`.ProgressPlotter`
        The plotter to test.

    path : :class:`pathlib2.Path`
        The path into which the plot is written.

    plot_data : :class:`pandas.DataFrame`
        The correct plot data.

    Returns
    -------
    None : :class:`NoneType`

    """

    plotter.write(path)
    assert plotter.get_plot_data().equals(plot_data)
