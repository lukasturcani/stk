"""
Dibormo
=======

"""

from __future__ import annotations

import typing

from .utilities import get_atom_map
from .generic_functional_group import GenericFunctionalGroup
from ...atoms import Br, Atom


__all__ = (
    'Dibromo',
)


class Dibromo(GenericFunctionalGroup):
    """
    Represents a dibromo functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[bromine1][atom1][atom2][bromine2]``.

    """

    def __init__(
        self,
        bromine1: Br,
        atom1: Atom,
        bromine2: Br,
        atom2: Atom,
        bonders: tuple[Atom, ...],
        deleters: tuple[Atom, ...],
        placers: typing.Optional[tuple[Atom, ...]] = None,
    ) -> None:
        """
        Initialize a :class:`.Dibromo` instance.

        Parameters:

            bromine1:
                The ``[bromine1]`` atom.

            atom1:
                The ``[atom1]`` atom.

            bromine2:
                The ``[bromine2]`` atom.

            atom2:
                The ``[atom]`` atom.

            bonders:
                The bonder atoms.

            deleters:
                The deleter atoms.

            placers:
                The placer atoms. If ``None`` the `bonders` will be
                used.

        """

        GenericFunctionalGroup.__init__(
            self=self,
            atoms=(bromine1, atom1, bromine2, atom2),
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )
        self._bromine1 = bromine1
        self._atom1 = atom1
        self._bromine2 = bromine2
        self._atom2 = atom2

    def get_atom1(self) -> Atom:
        """
        Get the ``[atom1]`` atom.

        Returns:

            The ``[atom1]`` atom.

        """

        return self._atom1

    def get_bromine1(self) -> Br:
        """
        Get the ``[bromine1]`` atom.

        Returns:

            The ``[bromine1]`` atom.

        """

        return self._bromine1

    def get_atom2(self) -> Atom:
        """
        Get the ``[atom2]`` atom.

        Returns:

            The ``[atom2]`` atom.

        """

        return self._atom2

    def get_bromine2(self) -> Br:
        """
        Get the ``[bromine2]`` atom.

        Returns:

            The ``[bromine2]`` atom.

        """

        return self._bromine2

    def clone(self) -> Dibromo:
        clone = self._clone()
        clone._atom1 = self._atom1
        clone._bromine1 = self._bromine1
        clone._atom2 = self._atom2
        clone._bromine2 = self._bromine2
        return clone

    def with_ids(
        self,
        id_map: dict[int, int],
    ) -> Dibromo:

        atom_map = get_atom_map(
            id_map=id_map,
            atoms=(
                *self._atoms,
                *self._placers,
                *self._core_atoms,
                *self._bonders,
                *self._deleters,
                self._atom1,
                self._bromine1,
                self._atom2,
                self._bromine2,
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
        clone._atom1 = atom_map.get(
            self._atom1.get_id(),
            self._atom1,
        )
        clone._bromine1 = atom_map.get(
            self._bromine1.get_id(),
            self._bromine1,
        )
        clone._atom2 = atom_map.get(
            self._atom2.get_id(),
            self._atom2,
        )
        clone._bromine2 = atom_map.get(
            self._bromine2.get_id(),
            self._bromine2,
        )
        return clone

    def __repr__(self) -> str:
        return (
            f'{self.__class__.__name__}('
            f'{self._bromine1}, {self._atom1}, {self._bromine2}, '
            f'{self._atom2}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
