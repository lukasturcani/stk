"""
Boronic Acid
============

"""

from __future__ import annotations

import typing

from . import utilities as _utilities
from . import generic_functional_group as _generic_functional_group
from ... import atoms as _atoms

__all__ = (
    'BoronicAcid',
)


class BoronicAcid(
    _generic_functional_group.GenericFunctionalGroup,
):
    """
    Represents a boronic acid functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][boron]([oxygen1][hydrogen1])[oxygen2][hydrogen2]``.

    """

    def __init__(
        self,
        boron: _atoms.B,
        # O is not an ambiguous name.
        oxygen1: _atoms.O,  # noqa
        hydrogen1: _atoms.H,
        # O is not an ambiguous name.
        oxygen2: _atoms.O,  # noqa
        hydrogen2: _atoms.H,
        atom: _atoms.Atom,
        bonders: tuple[_atoms.Atom, ...],
        deleters: tuple[_atoms.Atom, ...],
        placers: typing.Optional[tuple[_atoms.Atom, ...]] = None,
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
                the bonder atoms.

            deleters: :class:`tuple` of :class:`.Atom`
                The deleter atoms.

            placers:
                The placer atoms. If ``None`` the `bonders` will be
                used.

        """

        _generic_functional_group.GenericFunctionalGroup.__init__(
            self=self,
            atoms=(
                boron,
                oxygen1,
                hydrogen1,
                oxygen2,
                hydrogen2,
                atom,
            ),
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )
        self._boron = boron
        self._oxygen1 = oxygen1
        self._hydrogen1 = hydrogen1
        self._oxygen2 = oxygen2
        self._hydrogen2 = hydrogen2
        self._atom = atom

    def get_boron(self) -> _atoms.B:
        """
        Get the ``[boron]`` atom.

        Returns:

            The ``[boron]`` atom.

        """

        return self._boron

    # O is not an ambiguous name.
    def get_oxygen1(self) -> _atoms.O:  # noqa
        """
        Get the ``[oxygen1]`` atom.

        Returns:

            The ``[oxygen1]`` atom.

        """

        return self._oxygen1

    def get_hydrogen1(self) -> _atoms.H:
        """
        Get the ``[hydrogen1]`` atom.

        Returns:

            The ``[hydrogen1]`` atom.

        """

        return self._hydrogen1

    # O is not an ambiguous name.
    def get_oxygen2(self) -> _atoms.O:  # noqa
        """
        Get the ``[oxygen2]`` atom.

        Returns:

            The ``[oxygen2]`` atom.

        """

        return self._oxygen2

    def get_hydrogen2(self) -> _atoms.H:
        """
        Get the ``[hydrogen2]`` atom.

        Returns:

            The ``[hydrogen2]`` atom.

        """

        return self._hydrogen2

    def get_atom(self) -> _atoms.Atom:
        """
        Get the ``[atom]`` atom.

        Returns:

            The ``[atom]`` atom.

        """

        return self._atom

    def clone(self) -> BoronicAcid:
        clone = self._clone()
        clone._boron = self._boron
        clone._oxygen1 = self._oxygen1
        clone._hydrogen1 = self._hydrogen1
        clone._oxygen2 = self._oxygen2
        clone._hydrogen2 = self._hydrogen2
        clone._atom = self._atom
        return clone

    def with_ids(
        self,
        id_map: dict[int, int],
    ) -> BoronicAcid:

        atom_map = _utilities.get_atom_map(
            id_map=id_map,
            atoms=(
                *self._atoms,
                *self._placers,
                *self._core_atoms,
                *self._bonders,
                *self._deleters,
                self._boron,
                self._oxygen1,
                self._hydrogen1,
                self._oxygen2,
                self._hydrogen2,
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
        clone._boron = atom_map.get(
            self._boron.get_id(),
            self._boron,
        )
        clone._oxygen1 = atom_map.get(
            self._oxygen1.get_id(),
            self._oxygen1,
        )
        clone._hydrogen1 = atom_map.get(
            self._hydrogen1.get_id(),
            self._hydrogen1,
        )
        clone._oxygen2 = atom_map.get(
            self._oxygen2.get_id(),
            self._oxygen2,
        )
        clone._hydrogen2 = atom_map.get(
            self._hydrogen2.get_id(),
            self._hydrogen2,
        )
        clone._atom = atom_map.get(
            self._atom.get_id(),
            self._atom,
        )
        return clone

    def __repr__(self) -> str:
        return (
            f'{self.__class__.__name__}('
            f'{self._boron}, {self._oxygen1}, {self._hydrogen1}, '
            f'{self._oxygen2}, {self._hydrogen2}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
