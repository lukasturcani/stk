"""
Periodic Info
=============

Class holding periodic cell information.

"""

import logging
import numpy as np
from ...utilities import abs_cap

logger = logging.getLogger(__name__)


class PeriodicInfo:
    """
    Periodic cell information for periodic systems.

    """

    def __init__(self, cell_matrix):
        """
        Initialize a :class:`PeriodicInfo` instance.

        Parameters
        ----------
        cell_matrix : :class:`tuple` of :class:`np.array`
            Tuple of cell lattice vectors (shape: (3,)) in Angstrom.

        """

        self._cell_matrix = cell_matrix

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
            cell_matrix=self._cell_matrix
        )
        return clone

    def get_cell_matrix(self):
        return self._cell_matrix

    def _with_cell_matrix(self, cell_matrix):
        self._cell_matrix = cell_matrix
        return self

    def with_cell_matrix(self, cell_matrix):
        """
        Return a clone with cell_matrix.

        Parameters
        ----------
        cell_matrix : :class:`tuple` of :class:`np.array`
            Tuple of cell lattice vectors (shape: (3,)) in Angstrom.

        Returns
        -------
        :class:`.PeriodicInfo`
            The clone. Has the same cell as the original.

        """
        return self.clone()._with_cell_matrix(cell_matrix)

    def get_lengths_and_angles(self):
        """
        Convert cell matrix to lengths and angles.

        Lengths are in Angstrom. Angles are in degrees.

        This code is modified from the pymatgen source code [1]_.

        Returns
        -------
        :class:`tuple` of :class:`float`
            Lengths and angles that define cell matrix,

        References
        ----------
        .. [1] https://pymatgen.org/_modules/pymatgen/core/lattice.html

        """

        a, b, c = tuple(
            np.sqrt(np.sum(i ** 2)).tolist() for i in self._cell_matrix
        )

        m = self._cell_matrix
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

        return (a, b, c, alpha, beta, gamma)

    def __str__(self):
        a, b, c, alpha, beta, gamma = self.get_lengths_and_angles()
        return (
            f'{self.__class__.__name__}(a={a}, b={b}, '
            f'c={c}, alpha={alpha}, beta={beta}, '
            f'gamma={gamma})'
        )

    def __repr__(self):
        return str(self)
