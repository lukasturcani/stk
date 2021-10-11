"""
Fluoro
======

"""

from __future__ import annotations

import typing

from .utilities import get_atom_map
from .generic_functional_group import GenericFunctionalGroup
from ...atom import Atom
from ...elements import F


__all__ = (
    'Fluoro',
)


class Fluoro(GenericFunctionalGroup):
    """
    Represents a fluoro functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[fluorine][atom]``.

    """

    def __init__(
        self,
        fluorine: F,
        atom: Atom,
        bonders: tuple[Atom, ...],
        deleters: tuple[Atom, ...],
        placers: typing.Optional[tuple[Atom, ...]] = None,
    ) -> None:
        """
        Initialize a :class:`.Fluoro` instance.

        Parameters:

            fluorine:
                The ``[fluorine]`` atom.

            atom:
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
            atoms=(fluorine, atom),
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )
        self._fluorine = fluorine
        self._atom = atom

    def get_fluorine(self) -> F:
        """
        Get the ``[fluorine]`` atom.

        Returns:

            The ``[fluorine]`` atom.

        """

        return self._fluorine

    def get_atom(self) -> Atom:
        """
        Get the ``[atom]`` atom.

        Returns:

            The ``[atom]`` atom.

        """

        return self._atom

    def clone(self) -> Fluoro:
        clone = self._clone()
        clone._fluorine = self._fluorine
        clone._atom = self._atom
        return clone

    def with_ids(
        self,
        id_map: dict[int, int],
    ) -> Fluoro:
        atom_map = get_atom_map(
            id_map=id_map,
            atoms=(
                *self._atoms,
                *self._placers,
                *self._core_atoms,
                *self._bonders,
                *self._deleters,
                self._atom,
                self._fluorine,
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
        clone._fluorine = atom_map.get(
            self._fluorine.get_id(),
            self._fluorine,
        )
        clone._atom = atom_map.get(
            self._atom.get_id(),
            self._atom,
        )
        return clone

    def __repr__(self) -> str:
        return (
            f'{self.__class__.__name__}('
            f'{self._fluorine}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
