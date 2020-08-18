"""
Periodic Info
=============

Class holding periodic cell information.

"""

import logging
import numpy as np
from ..utilities import abs_cap

logger = logging.getLogger(__name__)


class PeriodicInfo:
    """
    Periodic cell information for periodic systems.

    """

    def __init__(self, x_vector, y_vector, z_vector):
        """
        Initialize a :class:`.PeriodicInfo` instance.

        Converts cell matrix to lengths and angles, where lengths are
        in Angstrom and angles are in degrees. This code is modified
        from the pymatgen source code [1]_.

        Parameters
        ----------
        x_vector : :class:`numpy.ndarray`
            Cell lattice vector of shape (3, ) in x direction in
            Angstrom.

        y_vector : :class:`numpy.ndarray`
            Cell lattice vector of shape (3, ) in y direction in
            Angstrom.

        z_vector : :class:`numpy.ndarray`
            Cell lattice vector of shape (3, ) in z direction in
            Angstrom.

        References
        ----------
        .. [1] https://pymatgen.org/_modules/pymatgen/core/lattice.html

        """
        self._x_vector = x_vector
        self._y_vector = y_vector
        self._z_vector = z_vector
        self._cell_matrix = (x_vector, y_vector, z_vector)

        a, b, c = tuple(
            np.sqrt(np.sum(i ** 2)).tolist() for i in self._cell_matrix
        )
        self._a = a
        self._b = b
        self._c = c

        lengths = (a, b, c)
        angles = np.zeros(3)
        for i in range(3):
            j = (i + 1) % 3
            k = (i + 2) % 3
            angles[i] = abs_cap(
                np.dot(self._cell_matrix[j], self._cell_matrix[k])
                / (lengths[j] * lengths[k])
            )
        angles = np.arccos(angles) * 180.0 / np.pi

        alpha, beta, gamma = angles
        self._alpha = alpha
        self._beta = beta
        self._gamma = gamma

    def clone(self):
        """
        Return a clone.

        Returns
        -------
        :class:`.PeriodicInfo`
            The clone. Has the same cell as the original.

        """

        clone = self.__class__.__new__(self.__class__)
        PeriodicInfo.__init__(
            self=clone,
            x_vector=self._x_vector,
            y_vector=self._y_vector,
            z_vector=self._z_vector,
        )
        return clone

    def get_x_vector(self):
        return self._x_vector

    def get_y_vector(self):
        return self._y_vector

    def get_z_vector(self):
        return self._z_vector

    def get_cell_matrix(self):
        return self._cell_matrix

    def get_a(self):
        return self._a

    def get_b(self):
        return self._b

    def get_c(self):
        return self._c

    def get_alpha(self):
        return self._alpha

    def get_beta(self):
        return self._beta

    def get_gamma(self):
        return self._gamma

    def __str__(self):

        return (
            f'{self.__class__.__name__}(a={self._a}, b={self._b}, '
            f'c={self._c}, alpha={self._alpha}, beta={self._beta}, '
            f'gamma={self._gamma})'
        )

    def __repr__(self):
        return str(self)
