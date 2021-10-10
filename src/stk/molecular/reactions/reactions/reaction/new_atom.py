"""
New Atom
========

"""

from __future__ import annotations

import numpy as np

from ....atoms import Atom


__all__ = (
    'NewAtom',
)


class NewAtom:
    """
    Represents a new atom added by a :class:`.Reaction`.

    """

    __slots__ = ['_atom', '_position']

    def __init__(
        self,
        atom: Atom,
        position: np.ndarray,
    ) -> None:
        """
        Initialize a :class:`.NewAtom` instance.

        Parameters:

            atom:
                The atom added by the reaction. Must have a negative
                id.

            position:
                The position of the new atom.

        """

        self._atom = atom
        self._position = np.array(position, dtype=np.float64)
        self._position.setflags(write=False)

    def get_atom(self) -> Atom:
        """
        Get the new atom.

        The new atom will have a negative id, and should be assigned
        a new one by the construction process.

        Returns:

            The new atom.

        """

        return self._atom

    def get_position(self) -> np.ndarray:
        """
        Get the position of the new atom.

        Returns:

            The position of the new atom.

        """

        return self._position
