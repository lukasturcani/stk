"""
New Atom
========

"""

import numpy as np


class NewAtom:
    """
    Represents a new atom added by a :class:`.Reaction`.

    """

    __slots__ = ['_atom', '_position']

    def __init__(self, atom, position):
        """
        Initialize a :class:`.NewAtom` instance.

        Parameters
        ----------
        atom : :class:`.Atom`
            The atom added by the reaction. Must have a negative id.

        position : :class:`numpy.ndarray`
            The position of the new atom.

        """

        self._atom = atom
        self._position = np.array(position, dtype=np.float64)
        self._position.setflags(write=False)

    def get_atom(self):
        """
        Get the new atom.

        The new atom will have a negative id, and should be assigned
        a new one by the construction process.

        Returns
        -------
        :class:`.Atom`
            The new atom.

        """

        return self._atom

    def get_position(self):
        """
        Get the position of the new atom.

        Returns
        -------
        :class:`numpy.ndarray`
            The position of the new atom.

        """

        return self._position
