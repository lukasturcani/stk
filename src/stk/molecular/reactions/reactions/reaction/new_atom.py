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
    *Unpacking the Atom and Position*

    The class can be elegantly unpacked into the atom and its position.

    .. testcode:: unpacking-the-atom-and-position

        import stk
        import numpy as np

        new_atoms = (
            stk.NewAtom(stk.C(-1), np.array([1., 2., 3.])),
            stk.NewAtom(stk.H(-2), np.array([2., 5., 7.])),
            stk.NewAtom(stk.F(-3), np.array([8., 2., 3.])),
        )
        for atom, position in new_atoms:
            # Do stuff.
            print(atom, position)

    .. testoutput:: unpacking-the-atom-and-position

        C(-1) [1. 2. 3.]
        H(-2) [2. 5. 7.]
        F(-3) [8. 2. 3.]

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
