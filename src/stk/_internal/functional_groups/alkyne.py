"""
Alkyne
======

"""

from __future__ import annotations

import typing

from stk._internal.atom import Atom
from stk._internal.elements import C

from .generic_functional_group import GenericFunctionalGroup


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
        placers: typing.Optional[tuple[Atom, ...]] = None,
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

        self._carbon1 = carbon1
        self._atom1 = atom1
        self._carbon2 = carbon2
        self._atom2 = atom2
        atoms = (carbon1, atom1, carbon2, atom2)
        super().__init__(
            atoms=atoms,
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )

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
        clone = super()._clone()
        clone._carbon1 = self._carbon1
        clone._atom1 = self._atom1
        clone._carbon2 = self._carbon2
        clone._atom2 = self._atom2
        return clone

    def _with_ids(
        self,
        id_map: dict[int, int],
    ) -> Alkyne:
        super()._with_ids(id_map)
        if (carbon1_id := self._carbon1.get_id()) in id_map:
            self._carbon1 = self._carbon1.with_id(id_map[carbon1_id])

        if (atom1_id := self._atom1.get_id()) in id_map:
            self._atom1 = self._atom1.with_id(id_map[atom1_id])

        if (carbon2_id := self._carbon2.get_id()) in id_map:
            self._carbon2 = self._carbon2.with_id(id_map[carbon2_id])

        if (atom2_id := self._atom2.get_id()) in id_map:
            self._atom2 = self._atom2.with_id(id_map[atom2_id])

        return self

    def with_ids(
        self,
        id_map: dict[int, int],
    ) -> Alkyne:
        return self.clone()._with_ids(id_map)

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"{self._carbon1}, {self._atom1}, {self._carbon2}, "
            f"{self._atom2}, "
            f"bonders={self._bonders}, deleters={self._deleters})"
        )
