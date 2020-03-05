class CaseData:
    """
    A test case.

    Attributes
    ----------
    molecule : :class:`.Molecule`
        The molecule to be tested.

    smiles : :class:`str`
        The canonical smiles for the molecule.

    """

    def __init__(self, molecule, smiles):
        self.molecule = molecule
        self.smiles = smiles
