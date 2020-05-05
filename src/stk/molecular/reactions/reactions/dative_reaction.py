"""
Dative Reaction
===============

"""

from itertools import chain

from .reaction import Reaction
from ...bonds import Bond


class NoMetalAtomError(Exception):
    ...


class DativeReaction(Reaction):
    """
    A reaction between two functional groups.

    The reaction creates a dative bond between the *bonder* atoms, and
    deletes any *deleter* atoms. Importantly, the direction of the bond
    is set such that `bonder2` is the metal atom.

    """

    def __init__(self, reaction):
        """
        Initialize a :class:`.DativeReaction` instance.

        Parameters
        ----------
        reaction : :class:`.Reaction`
            The reaction.

        """

        self._reaction = reaction

    def _get_new_atoms(self):
        return self._reaction._get_new_atoms()

    def _get_bond_directionality(self):
        """
        Get the correct bond direction for a dative bond.

        Dative bonds go: `organic` -> `metal`, where only the valency
        of `metal` is impacted by the bond.

        Raises
        ------
        :class:`NoMetalAtomError`
            If neither atom is a metal atom.

        """

        def is_metal_atom(atom):
            """
            Check if atom is a metal atom.

            Parameters
            ----------
            atom : :class:`.Atom`
                An atom.

            Returns
            -------
            :class:`bool`
                ``True`` if `atom` is a metal atom.

            """

            # Metal atomic numbers.
            metal_atomic_numbers = set(chain(
                list(range(21, 31)),
                list(range(39, 49)),
                list(range(72, 81))
            ))

            return atom.get_atomic_number() in metal_atomic_numbers

        bondera = next(self._reaction._functional_group1.get_bonders())
        bonderb = next(self._reaction._functional_group2.get_bonders())

        if is_metal_atom(bondera):
            bonder2 = bondera
            bonder1 = bonderb
        elif is_metal_atom(bonderb):
            bonder2 = bonderb
            bonder1 = bondera
        else:
            raise NoMetalAtomError(
                f'{bondera} and {bonderb} are metal atoms, so a dative'
                ' bond is not necessary.'
            )
        return bonder1, bonder2

    def _get_new_bonds(self):
        bonds = self._reaction._get_new_bonds()
        for bond in bonds:
            bonder1, bonder2 = self._get_bond_directionality()
            yield Bond(
                atom1=bonder1,
                atom2=bonder2,
                order=self._reaction._bond_order,
                periodicity=self._reaction._periodicity,
            )

    def _get_deleted_atoms(self):
        return self._reaction._get_deleted_atoms()
