import numpy as np
import stk


class CaseData:
    """
    A test case.

    Attributes:

        molecule:
            The molecule to be tested.

        smiles:
            The smiles for the molecule.

        name:
            The name of the test case.

        position_matrix:
            The position matrix of the molecule.

    """

    position_matrix: np.ndarray

    def __init__(
        self,
        molecule: stk.Molecule,
        smiles: str,
        name: str,
    ) -> None:
        self.molecule = molecule
        self.smiles = smiles
        self.name = name
        self.position_matrix = molecule.get_position_matrix()

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
        CaseData.__init__(clone, self.molecule, self.smiles, self.name)
        clone.position_matrix = np.load(path)
        return clone

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return (
            f"{self.__class__.__name__}("
            f"{self.molecule}, {self.smiles!r}, {self.name!r})"
        )
