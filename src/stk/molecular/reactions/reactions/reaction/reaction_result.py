"""
Reaction Result
===============

"""


class ReactionResult:
    """
    The result of a reaction.

    """

    __slots__ = [
        '_new_atoms',
        '_new_bonds',
        '_deleted_atoms',
        '_deleted_bonds',
    ]

    def __init__(
        self,
        new_atoms,
        new_bonds,
        deleted_atoms,
        deleted_bonds,
    ):
        """
        Initialize a :class:`.ReactionResult` instance.

        Parameters
        ----------
        new_atoms : :class:`tuple` of :class:`.NewAtom`
            The new atoms added by the reaction.

        new_bonds : :class:`tuple` of :class:`.Bond`
            The bonds added by the reaction.

        deleted_atoms : :class:`tuple` of :class:`.Atom`
            The atoms deleted by the reaction.

        deleted_bonds : :class:`tuple` of :class:`.Bond`
            The bonds deleted by the reaction.

        """

        self._new_atoms = new_atoms
        self._new_bonds = new_bonds
        self._deleted_atoms = deleted_atoms
        self._deleted_bonds = deleted_bonds

    def get_new_atoms(self):
        """
        Get the new atoms added by the reaction.

        Returns
        -------
        :class:`tuple` of :class:`.NewAtom`
            The new atoms added by the reaction.

        """

        return self._new_atoms

    def get_new_bonds(self):
        """
        Get the new bonds added by the reaction.

        Returns
        -------
        :class:`tuple` of :class:`.Bond`
            The new bonds added by the reaction.

        """
        return self._new_bonds

    def get_deleted_atoms(self):
        """
        Get the atoms deleted by the reaction.

        Returns
        -------
        :class:`tuple` of :class:`.Atom`
            The atoms deleted by the reaction.

        """

        return self._deleted_atoms

    def get_deleted_bonds(self):
        """
        Get the bonds deleted by the reaction.

        Returns
        -------
        :class:`tuple` of :class:`.Bond`
            The bonds deleted by the reaction.

        """

        return self._deleted_bonds
