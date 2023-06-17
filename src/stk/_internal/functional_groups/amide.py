"""
Amide
=====

"""

from __future__ import annotations

import typing

from stk._internal.atom import Atom
from stk._internal.elements import C, H, N, O

from .generic_functional_group import GenericFunctionalGroup


class Amide(GenericFunctionalGroup):
    """
    Represents an amide functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][carbon](=[oxygen])[nitrogen]([hydrogen1])[hydrogen2]``.

    """

    def __init__(
        self,
        carbon: C,
        oxygen: O,  # noqa: Not an ambiguous name.
        nitrogen: N,
        hydrogen1: H,
        hydrogen2: H,
        atom: Atom,
        bonders: tuple[Atom, ...],
        deleters: tuple[Atom, ...],
        placers: typing.Optional[tuple[Atom, ...]] = None,
    ) -> None:
        """
        Initialize a :class:`.Amide` instance.

        Parameters:

            carbon:
                The ``[carbon]`` atom.

            oxygen:
                The ``[oxygen]`` atom.

            nitrogen:
                The ``[nitrogen]`` atom.

            hydrogen1:
                The ``[hydrogen1]`` atom.

            hydrogen2:
                The ``[hydrogen2]`` atom.

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

        self._carbon = carbon
        self._oxygen = oxygen
        self._nitrogen = nitrogen
        self._hydrogen1 = hydrogen1
        self._hydrogen2 = hydrogen2
        self._atom = atom
        atoms = (carbon, oxygen, nitrogen, hydrogen1, hydrogen2, atom)
        super().__init__(
            atoms=atoms,
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )

    def get_carbon(self) -> C:
        """
        Get the ``[carbon]`` atom.

        Returns:

            The ``[carbon]`` atom.

        """

        return self._carbon

    def get_oxygen(self) -> O:  # noqa: Not an ambiguous name.
        """
        Get the ``[oxygen]`` atom.

        Returns:

            The ``[oxygen]`` atom.

        """

        return self._oxygen

    def get_nitrogen(self) -> N:
        """
        Get the ``[nitrogen]`` atom.

        Returns:

            The ``[nitrogen]`` atom.

        """

        return self._nitrogen

    def get_hydrogen1(self) -> H:
        """
        Get the ``[hydrogen1]`` atom.

        Returns:

            The ``[hydrogen1]`` atom.

        """

        return self._hydrogen1

    def get_hydrogen2(self) -> H:
        """
        Get the ``[hydrogen2]`` atom.

        Returns:

            The ``[hydrogen2]`` atom.

        """

        return self._hydrogen2

    def get_atom(self) -> Atom:
        """
        Get the ``[atom]`` atom.

        Returns:

            The ``[atom]`` atom.

        """

        return self._atom

    def clone(self) -> Amide:
        clone = super()._clone()
        clone._carbon = self._carbon
        clone._oxygen = self._oxygen
        clone._nitrogen = self._nitrogen
        clone._hydrogen1 = self._hydrogen1
        clone._hydrogen2 = self._hydrogen2
        clone._atom = self._atom
        return clone

    def _with_ids(
        self,
        id_map: dict[int, int],
    ) -> Amide:
        super()._with_ids(id_map)

        if (carbon_id := self._carbon.get_id()) in id_map:
            self._carbon = self._carbon.with_id(id_map[carbon_id])

        if (oxygen_id := self._oxygen.get_id()) in id_map:
            self._oxygen = self._oxygen.with_id(id_map[oxygen_id])

        if (nitrogen_id := self._nitrogen.get_id()) in id_map:
            self._nitrogen = self._nitrogen.with_id(
                id=id_map[nitrogen_id],
            )

        if (hydrogen1_id := self._hydrogen1.get_id()) in id_map:
            self._hydrogen1 = self._hydrogen1.with_id(
                id=id_map[hydrogen1_id],
            )

        if (hydrogen2_id := self._hydrogen2.get_id()) in id_map:
            self._hydrogen2 = self._hydrogen2.with_id(
                id=id_map[hydrogen2_id],
            )

        if (atom_id := self._atom.get_id()) in id_map:
            self._atom = self._atom.with_id(id_map[atom_id])

        return self

    def with_ids(
        self,
        id_map: dict[int, int],
    ) -> Amide:
        return self.clone()._with_ids(id_map)

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"{self._carbon}, {self._oxygen}, {self._nitrogen}, "
            f"{self._hydrogen1}, {self._hydrogen2}, {self._atom}, "
            f"bonders={self._bonders}, deleters={self._deleters})"
        )
