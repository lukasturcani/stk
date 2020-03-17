class CaseData:
    """
    A test case.

    Attributes
    ----------
    molecule : :class:`.Molecule`
        The molecule to be tested.

    smiles : :class:`str`
        The smiles for the molecule.

    position_matrix : :class:`numpy.ndarray`
        The position matrix of the molecule.

    """

    def __init__(self, molecule, smiles, position_matrix):
        self.molecule = molecule
        self.smiles = smiles
        self.position_matrix = position_matrix

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self.molecule}, {self.smiles!r})'
        )
