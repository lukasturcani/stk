"""
PDB Updating Utilities
======================

"""

import rdkit.Chem.AllChem as rdkit

from stk.utilities import remake


def _with_structure_from_pdb(self, path):
    """
    Change structure to match a ``.pdb`` file.

    Parameters
    ----------
    path : :class:`str`
        The full path of the ``.pdb`` file from which the structure
        should be updated.

    Returns
    -------
    :class:`.Molecule`
        The molecule.

    """

    molecule = remake(
        rdkit.MolFromPDBFile(
            molFileName=path,
            sanitize=False,
            removeHs=False,
        )
    )
    return self._with_position_matrix(
        position_matrix=molecule.GetConformer().GetPositions()
    )
