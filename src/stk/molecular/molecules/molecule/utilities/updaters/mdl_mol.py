"""
MDL Mol Updating Utilities
==========================

"""

import rdkit.Chem.AllChem as rdkit

from stk.utilities import remake


def _with_structure_from_mol(self, path):
    """
    Change structure to match a ``.mol`` file.

    Parameters
    ----------
    path : :class:`str`
        The full path of the ``.mol`` file from which the structure
        should be updated.

    Returns
    -------
    :class:`.Molecule`
        The molecule.

    """

    molecule = remake(
        rdkit.MolFromMolFile(
            molFileName=path,
            sanitize=False,
            removeHs=False,
        )
    )
    return self._with_position_matrix(
        position_matrix=molecule.GetConformer().GetPositions()
    )
