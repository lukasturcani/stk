class CaseData:
    """
    A test case.

    Attributes
    ----------
    plotter : :class:`.ProgressPlotter`
        The plotter to test.

    plot_data : :class:`pandas.DataFrame`
        The correct plot data.

    """

    def __init__(self, plotter, plot_data):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        plotter : :class:`.ProgressPlotter`
            The plotter to test.

        plot_data : :class:`pandas.DataFrame`
            The correct plot data.

        """

        self.plotter = plotter
        self.plot_data = plot_data
