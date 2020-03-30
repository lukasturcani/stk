class CaseData:
    """
    A test case.

    Attributes
    ----------
    mutator : :class:`.MoleculeMutator` or \
            :class:`.ConstructedMoleculeMutator`
        The mutator to test.

    record : :class:`.MoleculeRecord` or \
            :class:`.ConstructedMoleculeRecord`
        The molecule to mutate.

    mutation_record : :class:`.MoleculeMutationRecord` or \
            :class:`.ConstructedMoleculeMutationRecord`
        The correct mutation record.

    """

    def __init__(self, mutator, record, mutation_record):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        mutator : :class:`.MoleculeMutator` or \
                :class:`.ConstructedMoleculeMutator`
            The mutator to test.

        record : :class:`.MoleculeRecord` or \
                :class:`.ConstructedMoleculeRecord`
            The molecule to mutate.

        mutation_record : :class:`.MoleculeMutationRecord` or \
                :class:`.ConstructedMoleculeMutationRecord`
            The correct mutation record.

        """

        self._mutator = mutator
        self._record = record
        self._mutation_record = mutation_record
