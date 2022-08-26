import stk


def test_selection_plotter(tmp_path, case_data):
    """
    Test :class:`.SelectionPlotter`.

    Parameters
    ----------
    tmp_path : :class:`pathlib2.Path`
        A directory into which the plot is written.

    case_data : :class:`.CaseData`
        A test case. Holds the selector whose selections will be
        plotted.

    Returns
    -------
    None : :class:`NoneType`

    """

    _test_selection_plotter(
        selector=case_data.selector,
        population=case_data.population,
        filename=tmp_path / "selection",
    )


def _test_selection_plotter(selector, population, filename):
    """
    Test :class:`.SelectionPlotter`.

    Parameters
    ----------
    selector : :class:`.Selector`
        The selector whose selections will be plotted.

    population : :class:`tuple` of :class:`.MoleculeRecord`
        The population from which `selector` selects.

    filename : :class:`str`
        The filename for the plots.

    Returns
    -------
    None : :class:`NoneType`

    """

    stk.SelectionPlotter(filename, selector)
    tuple(selector.select(population))
