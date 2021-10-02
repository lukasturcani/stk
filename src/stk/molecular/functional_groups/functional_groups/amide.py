"""
Amide
=====

"""

from __future__ import annotations

import typing

from . import utilities as _utilities
from . import generic_functional_group as _generic_functional_group
from ... import atoms as _atoms


__all__ = (
    'Amide',
)


class Amide(
    _generic_functional_group.GenericFunctionalGroup,
):
    """
    Represents an amide functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][carbon](=[oxygen])[nitrogen]([hydrogen1])[hydrogen2]``.

    """

    def __init__(
        self,
        carbon: _atoms.C,
        # O is not an ambiguous name.
        oxygen: _atoms.O,  # noqa
        nitrogen: _atoms.N,
        hydrogen1: _atoms.H,
        hydrogen2: _atoms.H,
        atom: _atoms.Atom,
        bonders: tuple[_atoms.Atom, ...],
        deleters: tuple[_atoms.Atom, ...],
        placers: typing.Optional[tuple[_atoms.Atom, ...]] = None,
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

        _generic_functional_group.GenericFunctionalGroup.__init__(
            self=self,
            atoms=(
                carbon,
                oxygen,
                nitrogen,
                hydrogen1,
                hydrogen2,
                atom,
            ),
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )
        self._carbon = carbon
        self._oxygen = oxygen
        self._nitrogen = nitrogen
        self._hydrogen1 = hydrogen1
        self._hydrogen2 = hydrogen2
        self._atom = atom

    def get_carbon(self) -> _atoms.C:
        """
        Get the ``[carbon]`` atom.

        Returns:

            The ``[carbon]`` atom.

        """

        return self._carbon

    # O is not an ambiguous name.
    def get_oxygen(self) -> _atoms.O:  # noqa
        """
        Get the ``[oxygen]`` atom.

        Returns:

            The ``[oxygen]`` atom.

        """

        return self._oxygen

    def get_nitrogen(self) -> _atoms.N:
        """
        Get the ``[nitrogen]`` atom.

        Returns:

            The ``[nitrogen]`` atom.

        """

        return self._nitrogen

    def get_hydrogen1(self) -> _atoms.H:
        """
        Get the ``[hydrogen1]`` atom.

        Returns:

            The ``[hydrogen1]`` atom.

        """

        return self._hydrogen1

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

    def clone(self) -> Amide:
        clone = self._clone()
        clone._carbon = self._carbon
        clone._oxygen = self._oxygen
        clone._nitrogen = self._nitrogen
        clone._hydrogen1 = self._hydrogen1
        clone._hydrogen2 = self._hydrogen2
        clone._atom = self._atom
        return clone

    def with_ids(
        self,
        id_map: dict[int, int],
    ) -> Amide:

        atom_map = _utilities.get_atom_map(
            id_map=id_map,
            atoms=(
                *self._atoms,
                *self._placers,
                *self._core_atoms,
                *self._bonders,
                *self._deleters,
                self._carbon,
                self._oxygen,
                self._nitrogen,
                self._hydrogen1,
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
        clone._carbon = atom_map.get(
            self._carbon.get_id(),
            self._carbon,
        )
        clone._oxygen = atom_map.get(
            self._oxygen.get_id(),
            self._oxygen,
        )
        clone._nitrogen = atom_map.get(
            self._nitrogen.get_id(),
            self._nitrogen,
        )
        clone._hydrogen1 = atom_map.get(
            self._hydrogen1.get_id(),
            self._hydrogen1,
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
            f'{self._carbon}, {self._oxygen}, {self._nitrogen}, '
            f'{self._hydrogen1}, {self._hydrogen2}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
