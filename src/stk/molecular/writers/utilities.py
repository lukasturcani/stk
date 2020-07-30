"""
This module defines utilities for writers.

"""

import numpy as np


def abs_cap(val, max_abs_val=1):
    """
    Returns the value with its absolute value capped at max_abs_val.

    Particularly useful in passing values to trignometric functions
    where numerical errors may result in an argument > 1 being passed
    in.

    This code is modified from the pymatgen source code [1]_.

    References
    ----------
    .. [1] https://pymatgen.org/pymatgen.util.num.html

    """

    return max(min(val, max_abs_val), -max_abs_val)


def cell_matrix_to_lengths_angles(periodic_cell):
    """
    Convert cell matrix to lengths and angles.

    Lengths are in Angstrom. Angles are in degrees.

    This code is modified from the pymatgen source code [1]_.

    References
    ----------
    .. [1] https://pymatgen.org/_modules/pymatgen/core/lattice.html

    """

    a, b, c = tuple(
        np.sqrt(np.sum(i ** 2)).tolist() for i in periodic_cell
    )

    m = periodic_cell
    lengths = (a, b, c)
    angles = np.zeros(3)
    for i in range(3):
        j = (i + 1) % 3
        k = (i + 2) % 3
        angles[i] = abs_cap(
            np.dot(m[j], m[k]) / (lengths[j] * lengths[k])
        )
    angles = np.arccos(angles) * 180.0 / np.pi

    alpha, beta, gamma = angles

    return [a, b, c, alpha, beta, gamma]
