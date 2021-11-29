"""
Periodic Info
=============

Class holding periodic cell information.

"""

import logging

import numpy as np

from ..utilities import cap_absolute_value

logger = logging.getLogger(__name__)


class PeriodicInfo:
    """
    Periodic cell information for periodic systems.

    """

    def __init__(self, vector_1, vector_2, vector_3):
        """
        Initialize a :class:`.PeriodicInfo` instance.

        Converts cell matrix to lengths and angles, where lengths are
        in Angstrom and angles are in degrees. This code is modified
        from the pymatgen source code [1]_.

        Parameters
        ----------
        vector_1 : :class:`numpy.ndarray`
            First cell lattice vector of shape (3, ) in Angstrom.

        vector_2 : :class:`numpy.ndarray`
            Second cell lattice vector of shape (3, ) in Angstrom.

        vector_3 : :class:`numpy.ndarray`
            Third cell lattice vector of shape (3, ) in Angstrom.

        References
        ----------
        .. [1] https://pymatgen.org/_modules/pymatgen/core/lattice.html

        """

        self._vector_1 = np.array(vector_1)
        self._vector_2 = np.array(vector_2)
        self._vector_3 = np.array(vector_3)
        self._cell_matrix = (
            self._vector_1,
            self._vector_2,
            self._vector_3,
        )

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
            angles[i] = cap_absolute_value(
                value=(
                    np.dot(
                        self._cell_matrix[j], self._cell_matrix[k]
                    ) / (lengths[j] * lengths[k])
                ),
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
            vector_1=self._vector_1,
            vector_2=self._vector_2,
            vector_3=self._vector_3,
        )
        return clone

    def get_vector_1(self):
        """
        Get *x* vector.

        Returns
        -------
        :class:`numpy.ndarray`
            Cell lattice vector of shape (3, ) in *x* direction in
            Angstrom.

        """

        return np.array(self._vector_1)

    def get_vector_2(self):
        """
        Get *y* vector.

        Returns
        -------
        :class:`numpy.ndarray`
            Cell lattice vector of shape (3, ) in *y* direction in
            Angstrom.

        """

        return np.array(self._vector_2)

    def get_vector_3(self):
        """
        Get *z* vector.

        Returns
        -------
        :class:`numpy.ndarray`
            Cell lattice vector of shape (3, ) in *z* direction in
            Angstrom.

        """

        return np.array(self._vector_3)

    def get_cell_matrix(self):
        """
        Get cell matrix.

        Returns
        -------
        :class:`tuple` of :class:`numpy.ndarray`
            Tuple of length three containing *x*, *y* and *z* direction
            lattice vector of shape (3, ) in Angstrom.

        """

        return tuple(map(np.array, self._cell_matrix))

    def get_a(self):
        """
        Get *a* length.

        Returns
        -------
        :class:`float`
            Length of cell along *a* direction in Angstrom.

        """

        return self._a

    def get_b(self):
        """
        Get *b* length.

        Returns
        -------
        :class:`float`
            Length of cell along *b* direction in Angstrom.

        """

        return self._b

    def get_c(self):
        """
        Get *c* length.

        Returns
        -------
        :class:`float`
            Length of cell along *c* direction in Angstrom.

        """

        return self._c

    def get_alpha(self):
        """
        Get *alpha* angle.

        Returns
        -------
        :class:`float`
            *Alpha* angle of cell in degrees.

        """

        return self._alpha

    def get_beta(self):
        """
        Get *beta* angle.

        Returns
        -------
        :class:`float`
            *Beta* angle of cell in degrees.

        """

        return self._beta

    def get_gamma(self):
        """
        Get *gamma* angle.

        Returns
        -------
        :class:`float`
            *Gamma* angle of cell in degrees.

        """

        return self._gamma

    def __str__(self):

        return (
            f'{self.__class__.__name__}(a={self._a}, b={self._b}, '
            f'c={self._c}, alpha={self._alpha}, beta={self._beta}, '
            f'gamma={self._gamma})'
        )

    def __repr__(self):
        return str(self)
