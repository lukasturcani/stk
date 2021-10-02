"""
MAE Updating Utilities
======================

"""


import typing
import pathlib
import numpy as np
from stk.utilities import mol_from_mae_file


def get_position_matrix_from_mae(
    path: typing.Union[pathlib.Path, str],
) -> np.ndarray:
    """
    Get the position matrix from an ``.mae`` file.

    Parameters:

        path:
            The full path to the ``.mae`` file which holds the
            position matrix.

    Returns:

        The position matrix.

    """

    molecule = mol_from_mae_file(path)
    return molecule.GetConformer().GetPositions()
