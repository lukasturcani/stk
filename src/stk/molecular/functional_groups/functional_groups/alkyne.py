"""
Alkyne
======

"""

from __future__ import annotations

from typing import Optional

from .utilities import get_atom_map
from .generic_functional_group import GenericFunctionalGroup
from ...atoms import Atom, C


class Alkyne(GenericFunctionalGroup):
    """
    Represents an alkyne functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[carbon1]([atom1])#[carbon2][atom2]``.

    """

    def __init__(
        self,
        carbon1: C,
        atom1: Atom,
        carbon2: C,
        atom2: Atom,
        bonders: tuple[Atom, ...],
        deleters: tuple[Atom, ...],
        placers: Optional[tuple[Atom, ...]] = None,
    ) -> None:
        """
        Initialize a :class:`.Alkyne` instance.

        Parameters:

            carbon1:
                The ``[carbon1]`` atom.

            atom1:
                The ``[atom1]`` atom.

            carbon2:
                The ``[carbon2]`` atom.

            atom2:
                The ``[atom2]`` atom.

            bonders:
                    The bonder atoms.

            deleters:
                The deleter atoms.

            placers:
                The placer atoms. If ``None`` the `bonders` will be
                used.

        """

        atoms = (carbon1, atom1, carbon2, atom2)
        GenericFunctionalGroup.__init__(
            self=self,
            atoms=atoms,
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )
        self._carbon1 = carbon1
        self._atom1 = atom1
        self._carbon2 = carbon2
        self._atom2 = atom2

    def get_atom1(self) -> Atom:
        """
        Get the ``[atom1]`` atom.

        Returns:

            The ``[atom1]`` atom.

        """

        return self._atom1

    def get_carbon1(self) -> C:
        """
        Get the ``[carbon1]`` atom.

        Returns:

            The ``[carbon1]`` atom.

        """

        return self._carbon1

    def get_carbon2(self) -> C:
        """
        Get the ``[carbon2]`` atom.

        Returns:

            The ``[carbon2]`` atom.

        """

        return self._carbon2

    def get_atom2(self) -> Atom:
        """
        Get the ``[atom2]`` atom.

        Returns:

            The ``[atom2]`` atom.

        """

        return self._atom2

    def clone(self) -> Alkyne:
        clone = self._clone()
        clone._carbon1 = self._carbon1
        clone._atom1 = self._atom1
        clone._carbon2 = self._carbon2
        clone._atom2 = self._atom2
        return clone

    def with_ids(
        self,
        id_map: dict[int, int],
    ) -> Alkyne:

        atom_map = get_atom_map(
            id_map=id_map,
            atoms=(
                *self._atoms,
                *self._placers,
                *self._core_atoms,
                *self._bonders,
                *self._deleters,
                self._carbon1,
                self._atom1,
                self._carbon2,
                self._atom2,
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
        clone._carbon1 = atom_map.get(
            self._carbon1.get_id(),
            self._carbon1,
        )
        clone._atom1 = atom_map.get(
            self._atom1.get_id(),
            self._atom1,
        )
        clone._carbon2 = atom_map.get(
            self._carbon2.get_id(),
            self._carbon2,
        )
        clone._atom2 = atom_map.get(
            self._atom2.get_id(),
            self._atom2,
        )
        return clone

    def __repr__(self) -> str:
        return (
            f'{self.__class__.__name__}('
            f'{self._carbon1}, {self._atom1}, {self._carbon2}, '
            f'{self._atom2}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
