class CaseData:
    """
    A test case.

    Attributes
    ----------
    crosser : :class:`.MoleculeCrosser`
        The crosser to test.

    records : :class:`tuple` of :class:`.MoleculeRecord`
        The molecules to cross.

    crossover_records : :class:`tuple` of :class:`.CrossoverRecord`
        The correct offspring.

    """

    def __init__(self, crosser, records, crossover_records):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        crosser : :class:`.MoleculeCrosser`
            The crosser to test.

        records : :class:`tuple` of :class:`.MoleculeRecord`
            The molecules to cross.

        crossover_records : :class:`tuple` :class:`.CrossoverRecord`
            The correct offspring.

        """

        self.crosser = crosser
        self.records = records
        self.crossover_records = crossover_records
