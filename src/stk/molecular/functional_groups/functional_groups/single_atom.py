"""
Single Atom
===========

"""

from __future__ import annotations

from .utilities import get_atom_map
from .generic_functional_group import GenericFunctionalGroup
from ...atoms import Atom


__all__ = (
    'SingleAtom',
)


class SingleAtom(GenericFunctionalGroup):
    """
    Represents an abstract single atom functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom]``.

    """

    def __init__(self, atom: Atom) -> None:
        """
        Initialize a :class:`.SingleAtom` instance.

        Parameters:

            atom:
                Any :class:`.Atom` will work.

        """

        GenericFunctionalGroup.__init__(
            self=self,
            atoms=(atom, ),
            bonders=(atom, ),
            deleters=(),
            placers=(atom, )
        )
        self._atom = atom

    def get_atom(self) -> Atom:
        """
        Get the atom that defines the functional group.

        Returns:

            The atom.

        """

        return self._atom

    def with_ids(
        self,
        id_map: dict[int, int],
    ) -> SingleAtom:

        atom_map = get_atom_map(
            id_map=id_map,
            atoms=(
                *self._atoms,
                *self._placers,
                *self._core_atoms,
                *self._bonders,
                *self._deleters,
                self._atom,
            ),
        )
        clone = self.__class__.__new__(self.__class__)
        clone._atoms = tuple(
            atom_map.get(atom.get_id(), atom)
            for atom in self._atoms
        )
        clone._placers = tuple(
            atom_map.get(atom.get_id(), atom)
            for atom in self._placers
        )
        clone._core_atoms = tuple(
            atom_map.get(atom.get_id(), atom)
            for atom in self._core_atoms
        )
        clone._bonders = tuple(
            atom_map.get(atom.get_id(), atom)
            for atom in self._bonders
        )
        clone._deleters = tuple(
            atom_map.get(atom.get_id(), atom)
            for atom in self._deleters
        )
        clone._atom = atom_map.get(
            self._atom.get_id(),
            self._atom,
        )
        return clone

    def clone(self) -> SingleAtom:
        clone = self._clone()
        clone._atom = self._atom
        return clone

    def __repr__(self) -> str:
        return (
            f'{self.__class__.__name__}({self._atom})'
        )
