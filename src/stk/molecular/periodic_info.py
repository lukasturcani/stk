"""
Periodic Info
=============

Class holding periodic cell information.

"""

from __future__ import annotations

import logging

import numpy as np

from ..utilities import cap_absolute_value

logger = logging.getLogger(__name__)


class PeriodicInfo:
    """
    Periodic cell information for periodic systems.

    """

    def __init__(
        self,
        vector_1: np.ndarray,
        vector_2: np.ndarray,
        vector_3: np.ndarray,
    ) -> None:
        """
        Initialize a :class:`.PeriodicInfo` instance.

        Converts cell matrix to lengths and angles, where lengths are
        in Angstrom and angles are in degrees. This code is modified
        from the pymatgen source code [1]_.

        Parameters:

            vector_1:
                First cell lattice vector of shape (3, ) in Angstrom.

            vector_2:
                Second cell lattice vector of shape (3, ) in Angstrom.

            vector_3:
                Third cell lattice vector of shape (3, ) in Angstrom.

        References:

            .. [1] https://pymatgen.org/_modules/pymatgen/core/\
lattice.html

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

    def clone(self) -> PeriodicInfo:
        """
        Return a clone.

        Returns:

            The clone.

        """

        clone = self.__class__.__new__(self.__class__)
        PeriodicInfo.__init__(
            self=clone,
            vector_1=self._vector_1,
            vector_2=self._vector_2,
            vector_3=self._vector_3,
        )
        return clone

    def get_vector_1(self) -> np.ndarray:
        """
        Get *x* vector.

        Returns:

            Cell lattice vector of shape (3, ) in *x* direction in
            Angstrom.

        """

        return np.array(self._vector_1)

    def get_vector_2(self) -> np.ndarray:
        """
        Get *y* vector.

        Returns:

            Cell lattice vector of shape (3, ) in *y* direction in
            Angstrom.

        """

        return np.array(self._vector_2)

    def get_vector_3(self) -> np.ndarray:
        """
        Get *z* vector.

        Returns:

            Cell lattice vector of shape (3, ) in *z* direction in
            Angstrom.

        """

        return np.array(self._vector_3)

    def get_cell_matrix(
        self,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Get cell matrix.

        Returns:

            Tuple of length three containing *x*, *y* and *z* direction
            lattice vector of shape (3, ) in Angstrom.

        """

        a, b, c = self._cell_matrix
        return (np.array(a), np.array(b), np.array(c))

    def get_a(self) -> float:
        """
        Get *a* length.

        Returns:

            Length of cell along *a* direction in Angstrom.

        """

        return self._a

    def get_b(self) -> float:
        """
        Get *b* length.

        Returns:

            Length of cell along *b* direction in Angstrom.

        """

        return self._b

    def get_c(self) -> float:
        """
        Get *c* length.

        Returns:

            Length of cell along *c* direction in Angstrom.

        """

        return self._c

    def get_alpha(self) -> float:
        """
        Get *alpha* angle.

        Returns:

            *Alpha* angle of cell in degrees.

        """

        return self._alpha

    def get_beta(self) -> float:
        """
        Get *beta* angle.

        Returns:

            *Beta* angle of cell in degrees.

        """

        return self._beta

    def get_gamma(self) -> float:
        """
        Get *gamma* angle.

        Returns:

            *Gamma* angle of cell in degrees.

        """

        return self._gamma

    def __str__(self) -> str:

        return (
            f'{self.__class__.__name__}(a={self._a}, b={self._b}, '
            f'c={self._c}, alpha={self._alpha}, beta={self._beta}, '
            f'gamma={self._gamma})'
        )

    def __repr__(self) -> str:
        return str(self)
