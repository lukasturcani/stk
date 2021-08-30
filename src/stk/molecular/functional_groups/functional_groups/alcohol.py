"""
Alcohol
=======

"""

from __future__ import annotations

from typing import Optional, TypeVar

from .generic_functional_group import GenericFunctionalGroup
from ...atoms import Atom, O, H


_T = TypeVar('_T', bound='GenericFunctionalGroup')


class Alcohol(GenericFunctionalGroup):
    """
    Represents an alcohol functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][oxygen][hydrogen]``.

    """

    def __init__(
        self,
        # O is not an ambiguous name.
        oxygen: O,  # noqa
        hydrogen: H,
        atom: Atom,
        bonders: tuple[Atom, ...],
        deleters: tuple[Atom, ...],
        placers: Optional[tuple[Atom, ...]] = None,
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

        atoms = (oxygen, hydrogen, atom)
        super().__init__(
            atoms=atoms,
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )
        self._oxygen = oxygen
        self._hydrogen = hydrogen
        self._atom = atom

    # O is not an ambiguous name.
    def get_oxygen(self) -> O:  # noqa
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
        clone = super()._clone()
        clone._oxygen = self._oxygen
        clone._hydrogen = self._hydrogen
        clone._atom = self._atom
        return clone

    def with_ids(
        self,
        id_map: dict[int, int],
    ) -> Alcohol:

        clone = super()._with_ids(id_map)
        clone._oxygen = id_map.get(
            self._oxygen.get_id(),
            self._oxygen,
        )
        clone._hydrogen = id_map.get(
            self._hydrogen.get_id(),
            self._hydrogen,
        )
        clone._atom = id_map.get(
            self._atom.get_id(),
            self._atom,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._oxygen}, {self._hydrogen}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
