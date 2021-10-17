"""
Reaction Result
===============

"""


from ...atom import Atom
from ...bond import Bond
from .new_atom import NewAtom

__all__ = (
    'ReactionResult',
)


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
        new_atoms: tuple[NewAtom, ...],
        new_bonds: tuple[Bond, ...],
        deleted_atoms: tuple[Atom, ...],
        deleted_bonds: tuple[Bond, ...],
    ) -> None:
        """
        Initialize a :class:`.ReactionResult` instance.

        Parameters:

            new_atoms:
                The new atoms added by the reaction.

            new_bonds:
                The bonds added by the reaction.

            deleted_atoms:
                The atoms deleted by the reaction.

            deleted_bonds:
                The bonds deleted by the reaction.

        """

        self._new_atoms = new_atoms
        self._new_bonds = new_bonds
        self._deleted_atoms = deleted_atoms
        self._deleted_bonds = deleted_bonds

    def get_new_atoms(self) -> tuple[NewAtom, ...]:
        """
        Get the new atoms added by the reaction.

        Returns:

            The new atoms added by the reaction.

        """

        return self._new_atoms

    def get_new_bonds(self) -> tuple[Bond, ...]:
        """
        Get the new bonds added by the reaction.

        Returns:

            The new bonds added by the reaction.

        """
        return self._new_bonds

    def get_deleted_atoms(self) -> tuple[Atom, ...]:
        """
        Get the atoms deleted by the reaction.

        Returns:

            The atoms deleted by the reaction.

        """

        return self._deleted_atoms

    def get_deleted_bonds(self) -> tuple[Bond, ...]:
        """
        Get the bonds deleted by the reaction.

        Returns

            The bonds deleted by the reaction.

        """

        return self._deleted_bonds
