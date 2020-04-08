class CaseData:
    """
    A test case.

    Attributes
    ----------
    mutator : :class:`.MoleculeMutator`
        The mutator to test.

    record : :class:`.MoleculeRecord`
        The molecule to mutate.

    mutation_record : :class:`.MutationRecord`
        The correct mutation record.

    """

    def __init__(self, mutator, record, mutation_record):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        mutator : :class:`.MoleculeMutator`
            The mutator to test.

        record : :class:`.MoleculeRecord`
            The molecule to mutate.

        mutation_record : :class:`.MutationRecord`
            The correct mutation record.

        """

        self.mutator = mutator
        self.record = record
        self.mutation_record = mutation_record
