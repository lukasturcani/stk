"""
XYZ Updating Utilities
======================

"""

import numpy as np

import typing
import pathlib
from stk.utilities import periodic_table


def get_position_matrix_from_xyz(
    path: typing.Union[pathlib.Path, str],
) -> np.ndarray:
    """
    Get the position matrix from an ``.xyz`` file.

    Parameters:

        path:
            The full path to the ``.mol`` file which holds the
            position matrix.

    Returns:

        The position matrix.

    """

    with open(path, 'r') as f:
        _, _, *content = f.readlines()

    # Save all the coords in the file.
    new_coords = []
    for i, line in enumerate(content):
        element, *coords = line.split()
        # Handle XYZ files with capitilisation of element symbols.
        element = element.title()
        if element.isnumeric():
            element = periodic_table[int(element)]

        new_coords.append([float(i) for i in coords])

    return np.array(new_coords)
