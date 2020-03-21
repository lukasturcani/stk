"""
InChIKey
========

"""

import rdkit.Chem.AllChem as rdkit

from .molecule_key import MoleculeKey


class InchiKey(MoleculeKey):
    """
    """

    def __init__(self):
        """
        Initialize an :class:`InchiKey` instance.

        """

        super().__init__(
            name='InChIKey',
            key=self._get_inchi_key,
        )

    @staticmethod
    def _get_inchi_key(molecule):
        """
        Get the InChIKey of `molecule`.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule whose InChIKey is needed.

        Returns
        -------
        :class:`str`
            The InChIKey.

        """

        return rdkit.MolToInchiKey(molecule.to_rdkit_mol())
