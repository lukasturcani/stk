"""
Thioacid
========

"""

from __future__ import annotations

import typing

from .utilities import get_atom_map
from .generic_functional_group import GenericFunctionalGroup
from ...atom import Atom
from ...elements import C, O, S, H


__all__ = (
    'Thioacid',
)


class Thioacid(GenericFunctionalGroup):
    """
    Represents a thioacid functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][carbon](=[oxygen])[sulfur][hydrogen]``.

    """

    def __init__(
        self,
        carbon: C,
        # O is not an ambiguous name.
        oxygen: O,  # noqa
        sulfur: S,
        hydrogen: H,
        atom: Atom,
        bonders: tuple[Atom, ...],
        deleters: tuple[Atom, ...],
        placers: typing.Optional[tuple[Atom, ...]] = None,
    ) -> None:
        """
        Initialize a :class:`.Thioacid` functional group.

        Parameters:

            carbon:
                The ``[carbon]`` atom.

            oxygen:
                The ``[oxygen]`` atom.

            sulfur:
                The ``[sulfur]`` atom.

            hydrogen:
                The ``[hydrogen]`` atom.

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
            atoms=(carbon, oxygen, sulfur, hydrogen, atom),
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )
        self._carbon = carbon
        self._oxygen = oxygen
        self._sulfur = sulfur
        self._hydrogen = hydrogen
        self._atom = atom

    def get_carbon(self) -> C:
        """
        Get the ``[carbon]`` atom.

        Returns:

            The ``[carbon]`` atom.

        """

        return self._carbon

    # O is not an ambiguous name.
    def get_oxygen(self) -> O:  # noqa
        """
        Get the ``[oxygen]`` atom.

        Returns:

            The ``[oxygen]`` atom.

        """

        return self._oxygen

    def get_sulfur(self) -> S:
        """
        Get the ``[sulfur]`` atom.

        Returns
        -------
        :class:`.S`
            The ``[sulfur]`` atom.

        """

        return self._sulfur

    def get_hydrogen(self) -> H:
        """
        Get the ``[hydrogen]`` atom.

        Returns:

            Get the ``[hydrogen]`` atom.

        """

        return self._hydrogen

    def get_atom(self) -> Atom:
        """
        Get the ``[atom]`` atom.

        Returns:

            The ``[atom]`` atom.

        """

        return self._atom

    def clone(self) -> Thioacid:
        clone = self._clone()
        clone._carbon = self._carbon
        clone._oxygen = self._oxygen
        clone._sulfur = self._sulfur
        clone._hydrogen = self._hydrogen
        clone._atom = self._atom
        return clone

    def with_ids(
        self,
        id_map: dict[int, int],
    ) -> Thioacid:
        atom_map = get_atom_map(
            id_map=id_map,
            atoms=(
                *self._atoms,
                *self._placers,
                *self._core_atoms,
                *self._bonders,
                *self._deleters,
                self._carbon,
                self._oxygen,
                self._sulfur,
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
        clone._carbon = atom_map.get(
            self._carbon.get_id(),
            self._carbon,
        )
        clone._oxygen = atom_map.get(
            self._oxygen.get_id(),
            self._oxygen,
        )
        clone._sulfur = atom_map.get(
            self._sulfur.get_id(),
            self._sulfur,
        )
        clone._hydrogen = atom_map.get(
            self._hydrogen.get_id(),
            self._hydrogen,
        )
        clone._atom = atom_map.get(
            self._atom.get_id(),
            self._atom,
        )
        return clone

    def __repr__(self) -> str:
        return (
            f'{self.__class__.__name__}('
            f'{self._carbon}, {self._oxygen}, {self._sulfur}, '
            f'{self._hydrogen}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
