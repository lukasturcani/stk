"""
Dative Reaction
===============

"""

from ....bonds import Bond
from ..reaction import Reaction
from .utilities import is_metal


class DativeReaction(Reaction):
    """
    A reaction between two functional groups.

    The reaction creates a dative bond between the *bonder* atoms, and
    deletes any *deleter* atoms. Importantly, the direction of the bond
    is set to run from the non-metal atom to the metal atom.

    """

    def __init__(self, reaction):
        """
        Initialize a :class:`.DativeReaction` instance.

        Parameters
        ----------
        reaction : :class:`.Reaction`
            A reaction which should be made dative.

        """

        self._reaction = reaction

    def _get_new_atoms(self):
        return self._reaction._get_new_atoms()

    def _get_new_bonds(self):
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

    def _get_deleted_atoms(self):
        return self._reaction._get_deleted_atoms()

    def _get_deleted_bonds(self):
        return self._reaction._get_deleted_bonds()
