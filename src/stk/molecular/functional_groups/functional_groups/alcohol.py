"""
Alcohol
=======

"""

from __future__ import annotations

from typing import Optional

from .utilities import get_atom_map
from .generic_functional_group import GenericFunctionalGroup
from ...atoms import Atom, O, H


class Alcohol(GenericFunctionalGroup):
    """
    Represents an alcohol functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][oxygen][hydrogen]``.

    """

    def __init__(
        self,
        # O is not an ambiguous name.
        oxygen: O,  # noqa
        hydrogen: H,
        atom: Atom,
        bonders: tuple[Atom, ...],
        deleters: tuple[Atom, ...],
        placers: Optional[tuple[Atom, ...]] = None,
    ) -> None:
        """
        Initialize a :class:`.Alcohol` instance.

        Parameters:

            oxygen:
                The oxygen atom.

            hydrogen:
                The hydrogen atom.

            atom:
                The atom to which the alcohol is attached.

            bonders:
                The bonder atoms.

            deleters:
                The deleter atoms.

            placers:
                The placer atoms. If ``None`` the `bonders` will be
                used.

        """

        atoms = (oxygen, hydrogen, atom)
        GenericFunctionalGroup.__init__(
            self=self,
            atoms=atoms,
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )
        self._oxygen = oxygen
        self._hydrogen = hydrogen
        self._atom = atom

    # O is not an ambiguous name.
    def get_oxygen(self) -> O:  # noqa
        """
        Get the oxygen atom.

        Returns:

            The oxygen atom.

        """

        return self._oxygen

    def get_hydrogen(self) -> H:
        """
        Get the hydrogen atom.

        Returns:

            The hydrogen atom.

        """

        return self._hydrogen

    def get_atom(self) -> Atom:
        """
        Get the atom to which the functional group is attached.

        Returns:

            The atom to which the functional group is attached.

        """

        return self._atom

    def clone(self) -> Alcohol:
        clone = self._clone()
        clone._oxygen = self._oxygen
        clone._hydrogen = self._hydrogen
        clone._atom = self._atom
        return clone

    def with_ids(
        self,
        id_map: dict[int, int],
    ) -> Alcohol:

        atom_map = get_atom_map(
            id_map=id_map,
            atoms=(
                *self._atoms,
                *self._placers,
                *self._core_atoms,
                *self._bonders,
                *self._deleters,
                self._oxygen,
                self._hydrogen,
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
        clone._oxygen = atom_map.get(
            clone._oxygen.get_id(),
            clone._oxygen,
        )
        clone._hydrogen = atom_map.get(
            clone._hydrogen.get_id(),
            clone._hydrogen,
        )
        clone._atom = atom_map.get(
            clone._atom.get_id(),
            clone._atom,
        )
        return clone

    def __repr__(self) -> str:
        return (
            f'{self.__class__.__name__}('
            f'{self._oxygen}, {self._hydrogen}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
