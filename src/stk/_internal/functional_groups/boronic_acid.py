"""
Boronic Acid
============

"""

from __future__ import annotations

import typing

from stk._internal.atom import Atom
from stk._internal.elements import B, H, O

from .generic_functional_group import GenericFunctionalGroup


class BoronicAcid(GenericFunctionalGroup):
    """
    Represents a boronic acid functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][boron]([oxygen1][hydrogen1])[oxygen2][hydrogen2]``.

    """

    def __init__(
        self,
        boron: B,
        oxygen1: O,  # noqa: Not an ambiguous name.
        hydrogen1: H,
        oxygen2: O,  # noqa: Not an ambiguous name.
        hydrogen2: H,
        atom: Atom,
        bonders: tuple[Atom, ...],
        deleters: tuple[Atom, ...],
        placers: typing.Optional[tuple[Atom, ...]] = None,
    ) -> None:
        """
        Initialize a :class:`.BoronicAcid` instance.

        Parameters:

            boron:
                The ``[boron]`` atom.

            oxygen1:
                The ``[oxygen1]`` atom.

            hydrogen1:
                The ``[hydrogen]`` atom.

            oxygen2:
                The ``[oyxgen2]`` atom.

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

        self._boron = boron
        self._oxygen1 = oxygen1
        self._hydrogen1 = hydrogen1
        self._oxygen2 = oxygen2
        self._hydrogen2 = hydrogen2
        self._atom = atom
        atoms = (boron, oxygen1, hydrogen1, oxygen2, hydrogen2, atom)
        super().__init__(
            atoms=atoms,
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )

    def get_boron(self) -> B:
        """
        Get the ``[boron]`` atom.

        Returns:

            The ``[boron]`` atom.

        """

        return self._boron

    def get_oxygen1(self) -> O:  # noqa: Not an ambiguous name.
        """
        Get the ``[oxygen1]`` atom.

        Returns:

            The ``[oxygen1]`` atom.

        """

        return self._oxygen1

    def get_hydrogen1(self) -> H:
        """
        Get the ``[hydrogen1]`` atom.

        Returns:

            The ``[hydrogen1]`` atom.

        """

        return self._hydrogen1

    def get_oxygen2(self) -> O:  # noqa: Not an ambiguous name.
        """
        Get the ``[oxygen2]`` atom.

        Returns:

            The ``[oxygen2]`` atom.

        """

        return self._oxygen2

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

    def clone(self) -> BoronicAcid:
        clone = super()._clone()
        clone._boron = self._boron
        clone._oxygen1 = self._oxygen1
        clone._hydrogen1 = self._hydrogen1
        clone._oxygen2 = self._oxygen2
        clone._hydrogen2 = self._hydrogen2
        clone._atom = self._atom
        return clone

    def _with_ids(
        self,
        id_map: dict[int, int],
    ) -> BoronicAcid:
        super()._with_ids(id_map)

        if (boron_id := self._boron.get_id()) in id_map:
            self._boron = self._boron.with_id(id_map[boron_id])

        if (oxygen1_id := self._oxygen1.get_id()) in id_map:
            self._oxygen1 = self._oxygen1.with_id(id_map[oxygen1_id])

        if (hydrogen1_id := self._hydrogen1.get_id()) in id_map:
            self._hydrogen1 = self._hydrogen1.with_id(
                id=id_map[hydrogen1_id],
            )

        if (oxygen2_id := self._oxygen2.get_id()) in id_map:
            self._oxygen2 = self._oxygen2.with_id(id_map[oxygen2_id])

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
    ) -> BoronicAcid:
        return self.clone()._with_ids(id_map)

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"{self._boron}, {self._oxygen1}, {self._hydrogen1}, "
            f"{self._oxygen2}, {self._hydrogen2}, {self._atom}, "
            f"bonders={self._bonders}, deleters={self._deleters})"
        )
