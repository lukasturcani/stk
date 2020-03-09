"""
New Atom
========

"""

import numpy as np


class NewAtom:
    """
    Represents a new atom added by a :class:`.Reaction`.

    Examples
    --------
    The class can be elegantly unpacked into the atom and its position.

    .. code-block:: python

        import stk
        import numpy as np

        new_atoms = (
            NewAtom(stk.C(-1), np.array([1., 2., 3.])),
            NewAtom(stk.H(-2), np.array([2., 5., 7.])),
            NewAtom(stk.F(-3), np.array([8., 2., 3.])),
        )
        for atom, position in new_atoms:
            # Do stuff.

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

    def __iter__(self):
        """
        Iterate over the atom and position.

        """

        yield self._atom
        yield self._position
