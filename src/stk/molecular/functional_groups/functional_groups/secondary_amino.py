"""
Secondary Amino
===============

"""

from __future__ import annotations

import typing

from .utilities import get_atom_map
from .generic_functional_group import GenericFunctionalGroup
from ...atom import Atom
from ...elements import N, H


__all__ = (
    'SecondaryAmino',
)


class SecondaryAmino(GenericFunctionalGroup):
    """
    Represents a secondary amino functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom1][nitrogen]([hydrogen])[atom2]``.

    """

    def __init__(
        self,
        nitrogen: N,
        hydrogen: H,
        atom1: Atom,
        atom2: Atom,
        bonders: tuple[Atom, ...],
        deleters: tuple[Atom, ...],
        placers: typing.Optional[tuple[Atom, ...]] = None,
    ) -> None:
        """
        Initialize a :class:`.SecondaryAmine` instance.

        Parameters:

            nitrogen:
                The ``[nitrogen]`` atom

            hydrogen:
                The ``[hydrogen]`` atom.

            atom1:
                The ``[atom]`` atom.

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
            atoms=(nitrogen, hydrogen, atom1, atom2),
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )
        self._nitrogen = nitrogen
        self._hydrogen = hydrogen
        self._atom1 = atom1
        self._atom2 = atom2

    def get_nitrogen(self) -> N:
        """
        Get the ``[nitrogen]`` atom.

        Returns:

            The ``[nitrogen]`` atom.

        """

        return self._nitrogen

    def get_hydrogen(self) -> H:
        """
        Get the ``[hydrogen]`` atom.

        Returns:

            The ``[hydrogen]`` atom.

        """

        return self._hydrogen

    def get_atom1(self) -> Atom:
        """
        Get the ``[atom1]`` atom.

        Returns:

            The ``[atom1]`` atom.

        """

        return self._atom1

    def get_atom2(self) -> Atom:
        """
        Get the ``[atom2]`` atom.

        Returns:

            The ``[atom2]`` atom.

        """

        return self._atom2

    def clone(self) -> SecondaryAmino:
        clone = self._clone()
        clone._nitrogen = self._nitrogen
        clone._hydrogen = self._hydrogen
        clone._atom1 = self._atom1
        clone._atom2 = self._atom2
        return clone

    def with_ids(
        self,
        id_map: dict[int, int],
    ):
        atom_map = get_atom_map(
            id_map=id_map,
            atoms=(
                *self._atoms,
                *self._placers,
                *self._core_atoms,
                *self._bonders,
                *self._deleters,
                self._nitrogen,
                self._hydrogen,
                self._atom1,
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
        clone._nitrogen = atom_map.get(
            self._nitrogen.get_id(),
            self._nitrogen,
        )
        clone._hydrogen = atom_map.get(
            self._hydrogen.get_id(),
            self._hydrogen,
        )
        clone._atom1 = atom_map.get(
            self._atom1.get_id(),
            self._atom1,
        )
        clone._atom2 = atom_map.get(
            self._atom2.get_id(),
            self._atom2,
        )
        return clone

    def __repr__(self) -> str:
        return (
            f'{self.__class__.__name__}('
            f'{self._nitrogen}, {self._hydrogen}, {self._atom1}, '
            f'{self._atom2}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
