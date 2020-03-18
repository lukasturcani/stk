class CaseData:
    """
    A test case.

    Attributes
    ----------
    molecule : :class:`.Molecule`
        The molecule to test.

    atoms : :class:`tuple`
        The correct atoms of the molecule.

    """

    def __init__(self, molecule, atoms):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule to test.

        atoms : :class:`.Atom`
            The correct atoms of the molecule.

        """

        self.molecule = molecule
        self.atoms = atoms
