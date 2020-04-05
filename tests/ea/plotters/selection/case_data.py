class CaseData:
    """
    A test case.

    Attributes
    ----------
    selector : :class:`.Selector`
        The selector whose selections will be plotted.

    population : :class:`tuple` of :class:`.MoleculeRecord`
        The population from which the :attr:`.selector` selects.


    """

    def __init__(self, selector, population):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        selector : :class:`.Selector`
            The selector whose selections will be plotted.

        population : :class:`tuple` of :class:`.MoleculeRecord`
            The population from which the `selector` selects.

        """

        self.selector = selector
        self.population = population
