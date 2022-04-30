"""
Alcohol
=======

"""

from __future__ import annotations

import typing

from ...atoms import Atom, H, O
from .generic_functional_group import GenericFunctionalGroup


class Alcohol(GenericFunctionalGroup):
    """
    Represents an alcohol functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][oxygen][hydrogen]``.

    """

    def __init__(
        self,
        oxygen: O,  # noqa: Not an ambiguous name.
        hydrogen: H,
        atom: Atom,
        bonders: tuple[Atom, ...],
        deleters: tuple[Atom, ...],
        placers: typing.Optional[tuple[Atom, ...]] = None,
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

        self._oxygen = oxygen
        self._hydrogen = hydrogen
        self._atom = atom
        atoms = (oxygen, hydrogen, atom)
        super().__init__(
            atoms=atoms,
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )

    def get_oxygen(self) -> O:  # noqa: Not an ambiguous name.
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

    def _with_ids(
        self,
        id_map: dict[int, int],
    ) -> Alcohol:

        clone = super()._with_ids(id_map)

        if (oxygen_id := self._oxygen.get_id()) in id_map:
            clone._oxygen = self._oxygen.with_id(id_map[oxygen_id])

        if (hydrogen_id := self._hydrogen.get_id()) in id_map:
            clone._hydrogen = self._hydrogen.with_id(
                id=id_map[hydrogen_id],
            )

        if (atom_id := self._atom.get_id()) in id_map:
            clone._atom = self._atom.with_id(id_map[atom_id])

        return clone

    def with_ids(
        self,
        id_map: dict[int, int],
    ) -> Alcohol:

        return self.clone()._with_ids(id_map)

    def __repr__(self) -> str:
        return (
            f'{self.__class__.__name__}('
            f'{self._oxygen}, {self._hydrogen}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
