"""
MDL Mol Updating Utilities
==========================

"""

import typing
import pathlib
import numpy as np
import rdkit.Chem.AllChem as rdkit

from stk.utilities import remake


def get_position_matrix_from_mol(
    path: typing.Union[pathlib.Path, str],
) -> np.ndarray:
    """
    Get the position matrix from a ``.mol`` file.

    Parameters:

        path:
            The full path to the ``.mol`` file which holds the
            position matrix.

    Returns:

        The position matrix.

    """

    molecule = remake(
        rdkit.MolFromMolFile(
            molFileName=path,
            sanitize=False,
            removeHs=False,
        )
    )
    return molecule.GetConformer().GetPositions()
