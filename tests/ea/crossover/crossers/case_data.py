class CaseData:
    """
    A test case.

    Attributes
    ----------
    crosser : :class:`.MoleculeCrosser` or \
            :class:`.ConstructedMoleculeCrosser`
        The crosser to test.

    records : :class:`tuple` of :class:`.MoleculeRecord` or \
            :class:`tuple` of :class:`.ConstructedMoleculeRecord`
        The molecules to cross.

    crossover_records : :class:`tuple` of \
            :class:`MoleculeCrossoverRecord` or \
            :class:`tuple` of \
            :class:`.ConstructedMoleculeCrossoverRecord`
        The correct offspring.

    """

    def __init__(self, crosser, records, crossover_records):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        crosser : :class:`.MoleculeCrosser` or \
                :class:`.ConstructedMoleculeCrosser`
            The crosser to test.

        records : :class:`tuple` of :class:`.MoleculeRecord` or \
                :class:`tuple` of :class:`.ConstructedMoleculeRecord`
            The molecules to cross.

        crossover_records : :class:`tuple` of \
                :class:`MoleculeCrossoverRecord` or \
                :class:`tuple` of \
                :class:`.ConstructedMoleculeCrossoverRecord`
            The correct offspring.

        """

        self.crosser = crosser
        self.records = records
        self.crossover_records = crossover_records
