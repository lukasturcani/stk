"""
Diol
====

"""

from __future__ import annotations

import typing

from .utilities import get_atom_map
from .generic_functional_group import GenericFunctionalGroup
from ...atoms import O, H, Atom


__all__ = (
    'Diol',
)


class Diol(GenericFunctionalGroup):
    """
    Represents a diol functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[hydrogen1][oxygen1][atom1][atom2][oxygen2][hydrogen2]``.

    """

    def __init__(
        self,
        atom1: Atom,
        # O is not an ambiguous name.
        oxygen1: O,  # noqa
        hydrogen1: H,
        atom2: Atom,
        # O is not an ambiguous name.
        oxygen2: O,  # noqa
        hydrogen2: H,
        bonders: tuple[Atom, ...],
        deleters: tuple[Atom, ...],
        placers: typing.Optional[tuple[Atom, ...]] = None,
    ):
        """
        Initialize a :class:`.Diol` instance.

        Parameters:

            atom1:
                The ``[atom1]`` atom.

            oxygen1:
                The ``[oxygen1]`` atom.

            hydrogen1:
                The ``[hydrogen1]`` atom.

            atom2:
                The ``[atom2]`` atom.

            oxygen2:
                The ``[oxygen2]`` atom.

            hydrogen2:
                The ``[hydrogen2]`` atom.

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
            atoms=(
                atom1,
                oxygen1,
                hydrogen1,
                atom2,
                oxygen2,
                hydrogen2,
            ),
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )
        self._atom1 = atom1
        self._oxygen1 = oxygen1
        self._hydrogen1 = hydrogen1
        self._atom2 = atom2
        self._oxygen2 = oxygen2
        self._hydrogen2 = hydrogen2

    def get_atom1(self) -> Atom:
        """
        Get the ``[atom1]`` atom.

        Returns:

            The ``[atom1]`` atom.

        """

        return self._atom1

    # O is not an ambiguous name.
    def get_oxygen1(self) -> O:  # noqa
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

    def get_atom2(self) -> Atom:
        """
        Get the ``[atom2]`` atom.

        Returns:

            The ``[atom2]`` atom.

        """

        return self._atom2

    # O is not an ambiguous name.
    def get_oxygen2(self) -> O:  # noqa
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

    def clone(self) -> Diol:
        clone = self._clone()
        clone._atom1 = self._atom1
        clone._oxygen1 = self._oxygen1
        clone._hydrogen1 = self._hydrogen1
        clone._atom2 = self._atom2
        clone._oxygen2 = self._oxygen2
        clone._hydrogen2 = self._hydrogen2
        return clone

    def with_ids(
        self,
        id_map: dict[int, int],
    ) -> Diol:

        atom_map = get_atom_map(
            id_map=id_map,
            atoms=(
                *self._atoms,
                *self._placers,
                *self._core_atoms,
                *self._bonders,
                *self._deleters,
                self._atom1,
                self._oxygen1,
                self._atom2,
                self._oxygen2,
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
        clone._atom1 = atom_map.get(
            self._atom1.get_id(),
            self._atom1,
        )
        clone._oxygen1 = atom_map.get(
            self._oxygen1.get_id(),
            self._oxygen1,
        )
        clone._hydrogen1 = atom_map.get(
            self._hydrogen1.get_id(),
            self._hydrogen1,
        )
        clone._atom2 = atom_map.get(
            self._atom2.get_id(),
            self._atom2,
        )
        clone._oxygen2 = atom_map.get(
            self._oxygen2.get_id(),
            self._oxygen2,
        )
        clone._hydrogen2 = atom_map.get(
            self._hydrogen2.get_id(),
            self._hydrogen2,
        )
        return clone

    def __repr__(self) -> str:
        return (
            f'{self.__class__.__name__}('
            f'{self._atom1}, {self._oxygen1}, {self._hydrogen1}, '
            f'{self._atom2}, {self._oxygen2}, {self._hydrogen2}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
