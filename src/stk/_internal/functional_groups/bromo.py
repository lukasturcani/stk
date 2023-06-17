"""
Bromo
=====

"""

from __future__ import annotations

import typing

from stk._internal.atom import Atom
from stk._internal.elements import Br

from .generic_functional_group import GenericFunctionalGroup


class Bromo(GenericFunctionalGroup):
    """
    Represents a bromo functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[bromine][atom]``.

    """

    def __init__(
        self,
        bromine: Br,
        atom: Atom,
        bonders: tuple[Atom, ...],
        deleters: tuple[Atom, ...],
        placers: typing.Optional[tuple[Atom, ...]] = None,
    ) -> None:
        """
        Initialize a :class:`.Bromo` instance.

        Parameters:

            bromine:
                The ``[bromine]`` atom.

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

        self._bromine = bromine
        self._atom = atom
        super().__init__(
            atoms=(bromine, atom),
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )

    def get_bromine(self) -> Br:
        """
        Get the ``[bromine]`` atom.

        Returns:

            The ``[bromine]`` atom.

        """

        return self._bromine

    def get_atom(self) -> Atom:
        """
        Get the ``[atom]`` atom.

        Returns:

            The ``[atom]`` atom.

        """

        return self._atom

    def clone(self) -> Bromo:
        clone = super()._clone()
        clone._bromine = self._bromine
        clone._atom = self._atom
        return clone

    def _with_ids(
        self,
        id_map: dict[int, int],
    ) -> Bromo:
        super()._with_ids(id_map)
        if (bromine_id := self._bromine.get_id()) in id_map:
            self._bromine = self._bromine.with_id(id_map[bromine_id])

        if (atom_id := self._atom.get_id()) in id_map:
            self._atom = self._atom.with_id(id_map[atom_id])

        return self

    def with_ids(
        self,
        id_map: dict[int, int],
    ) -> Bromo:
        return self.clone()._with_ids(id_map)

    def __repr__(self) -> str:
        return (
            f"{self.__class__.__name__}("
            f"{self._bromine}, {self._atom}, "
            f"bonders={self._bonders}, deleters={self._deleters})"
        )
