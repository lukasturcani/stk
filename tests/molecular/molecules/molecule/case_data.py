import numpy as np


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

    def __init__(self, molecule, smiles):
        self.molecule = molecule
        self.smiles = smiles
        self.position_matrix = None

    def with_position_matrix(self, path):
        """
        Return a clone holding a position matrix.

        Parameters
        ----------
        path : :class:`str`
            The for to a dumped :class:`numpy.ndarray` position
            matrix.

        Returns
        -------
        :class:`.CaseData`
            The clone, holding the position matrix found in `path`.

        """

        clone = self.__class__.__new__(self.__class__)
        CaseData.__init__(clone, self.molecule, self.smiles)
        clone.position_matrix = np.load(path)
        return clone

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self.molecule}, {self.smiles!r})'
        )
