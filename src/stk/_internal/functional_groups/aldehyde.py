"""
Aldehyde
========

"""

from __future__ import annotations

import typing

from stk._internal.atom import Atom
from stk._internal.elements import C, H, O

from .generic_functional_group import GenericFunctionalGroup


class Aldehyde(GenericFunctionalGroup):
    """
    Represents an aldehyde functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][carbon](=[oxygen])[hydrogen]``.

    """

    def __init__(
        self,
        carbon: C,
        oxygen: O,  # noqa: Not an ambiguous name.
        hydrogen: H,
        atom: Atom,
        bonders: tuple[Atom, ...],
        deleters: tuple[Atom, ...],
        placers: typing.Optional[tuple[Atom, ...]] = None,
    ) -> None:
        """
        Initialize a :class:`.Aldehyde` instance.

        Parameters:

            carbon:
                The carbon atom.

            oxygen:
                The oxygen atom.

            hydrogen:
                The hydrogen atom.

            atom:
                The atom to which the functional group is attached.

            bonders:
                The bonder atoms.

            deleters:
                The deleter atoms.

            placers:
                The placer atoms. If ``None`` the `bonders` will be
                used.

        """

        self._carbon = carbon
        self._oxygen = oxygen
        self._hydrogen = hydrogen
        self._atom = atom
        atoms = (carbon, oxygen, hydrogen, atom)
        super().__init__(
            atoms=atoms,
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )

    def get_carbon(self) -> C:
        """
        Get the carbon atom.

        Returns:

            The carbon atom.

        """

        return self._carbon

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

    def _with_ids(
        self,
        id_map: dict[int, int],
    ) -> Aldehyde:
        super()._with_ids(id_map)

        if (carbon_id := self._carbon.get_id()) in id_map:
            self._carbon = self._carbon.with_id(id_map[carbon_id])

        if (oxygen_id := self._oxygen.get_id()) in id_map:
            self._oxygen = self._oxygen.with_id(id_map[oxygen_id])

        if (hydrogen_id := self._hydrogen.get_id()) in id_map:
            self._hydrogen = self._hydrogen.with_id(
                id=id_map[hydrogen_id],
            )

        if (atom_id := self._atom.get_id()) in id_map:
            self._atom = self._atom.with_id(id_map[atom_id])

        return self

    def with_ids(
        self,
        id_map: dict[int, int],
    ) -> Aldehyde:
        return self.clone()._with_ids(id_map)

    def clone(self) -> Aldehyde:
        clone = super()._clone()
        clone._carbon = self._carbon
        clone._oxygen = self._oxygen
        clone._hydrogen = self._hydrogen
        clone._atom = self._atom
        return clone

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"{self._carbon}, {self._oxygen}, {self._hydrogen}, "
            f"{self._atom}, bonders={self._bonders}, "
            f"deleters={self._deleters})"
        )
