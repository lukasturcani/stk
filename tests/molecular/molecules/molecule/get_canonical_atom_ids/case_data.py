class CaseData:
    """
    A test case.

    Attributes
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    canonical_atom_ids : :class:`dict`
        The correct mapping of atom ids in :attr:`.molecule` to their
        canonical atom ids.

    """

    def __init__(self, molecule, canonical_atom_ids):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule to test.

        canonical_atom_ids : :class:`dict`
            The correct mapping of atom ids in :attr:`.molecule` to
            their canonical atom ids.

        """

        self.molecule = molecule
        self.canonical_atom_ids = canonical_atom_ids
