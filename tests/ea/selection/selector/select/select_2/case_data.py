class CaseData:
    """
    A test case.

    Attributes
    ----------
    selector : :class:`.Selector`
        The selector to test.

    population : :class:`tuple` of :class:`.MoleculeRecord`
        The population from which batches are selected.

    selected : :class:`tuple` of :class:`.Batch`
            The batches which should be selected.

    """

    def __init__(self, selector, population, selected):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        selector : :class:`.Selector`
            The selector to test.

        population : :class:`tuple` of :class:`.MoleculeRecord`
            The population from which batches are selected.

        selected : :class:`tuple` of :class:`.Batch`
            The batches which should be selected.

        """

        self.selector = selector
        self.population = population
        self.selected = selected
