"""
Difluoro
========

"""

from __future__ import annotations

from typing import Optional

from .utilities import get_atom_map
from .generic_functional_group import GenericFunctionalGroup
from ...atoms import F, Atom


class Difluoro(GenericFunctionalGroup):
    """
    Represents a difluoro functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[fluorine1][atom1][atom2][fluorine2]``.

    """

    def __init__(
        self,
        fluorine1: F,
        atom1: Atom,
        fluorine2: F,
        atom2: Atom,
        bonders: tuple[Atom, ...],
        deleters: tuple[Atom, ...],
        placers: Optional[tuple[Atom, ...]] = None,
    ) -> None:
        """
        Initialize a :class:`.Difluoro` instance.

        Parameters:

            fluorine1:
                The ``[fluorine1]`` atom.

            atom1:
                The ``[atom1]`` atom.

            fluorine2:
                The ``[fluorine2]`` atom.

            atom2:
                The ``[atom2]`` atom.

            placers:
                The placer atoms. If ``None`` the `bonders` will be
                used.

        """

        GenericFunctionalGroup.__init__(
            self=self,
            atoms=(fluorine1, atom1, fluorine2, atom2),
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )
        self._fluorine1 = fluorine1
        self._atom1 = atom1
        self._fluorine2 = fluorine2
        self._atom2 = atom2

    def get_atom1(self) -> Atom:
        """
        Get the ``[atom1]`` atom.

        Returns:

            The ``[atom1]`` atom.

        """

        return self._atom1

    def get_fluorine1(self) -> F:
        """
        Get the ``[fluorine1]`` atom.

        Returns:

            The ``[fluorine1]`` atom.

        """

        return self._fluorine1

    def get_atom2(self) -> Atom:
        """
        Get the ``[atom2]`` atom.

        Returns:

            The ``[atom2]`` atom.

        """

        return self._atom2

    def get_fluorine2(self) -> F:
        """
        Get the ``[fluorine2]`` atom.

        Returns:

            The ``[fluorine2]`` atom.

        """

        return self._fluorine2

    def with_ids(
        self,
        id_map: dict[int, int],
    ) -> Difluoro:

        atom_map = get_atom_map(
            id_map=id_map,
            atoms=(
                *self._atoms,
                *self._placers,
                *self._core_atoms,
                *self._bonders,
                *self._deleters,
                self._atom1,
                self._fluorine1,
                self._atom2,
                self._fluorine2,
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
        clone._fluorine1 = atom_map.get(
            self._fluorine1.get_id(),
            self._fluorine1,
        )
        clone._atom2 = atom_map.get(
            self._atom2.get_id(),
            self._atom2,
        )
        clone._fluorine2 = atom_map.get(
            self._fluorine2.get_id(),
            self._fluorine2,
        )
        return clone

    def clone(self) -> Difluoro:
        clone = self._clone()
        clone._atom1 = self._atom1
        clone._fluorine1 = self._fluorine1
        clone._atom2 = self._atom2
        clone._fluorine2 = self._fluorine2
        return clone

    def __repr__(self) -> str:
        return (
            f'{self.__class__.__name__}('
            f'{self._fluorine1}, {self._atom1}, {self._fluorine2}, '
            f'{self._atom2}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
