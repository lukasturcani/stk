"""
Dative Reaction
===============

"""

from collections import abc

from .utilities import is_metal
from ..reaction import Reaction, NewAtom
from ....bonds import Bond
from ....atoms import Atom


__all__ = (
    'DativeReaction',
)


class DativeReaction(Reaction):
    """
    A reaction between two functional groups.

    The reaction creates a dative bond between the *bonder* atoms, and
    deletes any *deleter* atoms. Importantly, the direction of the bond
    is set to run from the non-metal atom to the metal atom.

    """

    def __init__(self, reaction: Reaction) -> None:
        """
        Initialize a :class:`.DativeReaction` instance.

        Parameters:

            reaction:
                A reaction which should be made dative.

        """

        self._reaction = reaction

    def _get_new_atoms(self) -> abc.Iterator[NewAtom]:
        return self._reaction._get_new_atoms()

    def _get_new_bonds(self) -> abc.Iterator[Bond]:
        for bond in self._reaction._get_new_bonds():
            if bond.get_order() == 9 and is_metal(bond.get_atom1()):
                yield Bond(
                    atom1=bond.get_atom2(),
                    atom2=bond.get_atom1(),
                    order=bond.get_order(),
                    periodicity=bond.get_periodicity(),
                )
            else:
                yield bond

    def _get_deleted_atoms(self) -> abc.Iterator[Atom]:
        return self._reaction._get_deleted_atoms()

    def _get_deleted_bonds(self) -> abc.Iterator[Bond]:
        return self._reaction._get_deleted_bonds()
