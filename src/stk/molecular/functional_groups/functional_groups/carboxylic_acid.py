"""
Carboxylic Acid
===============

"""

from __future__ import annotations

import typing

from . import utilities as _utilities
from . import generic_functional_group as _generic_functional_group
from ... import atoms as _atoms


__all__ = (
    'CarboxylicAcid',
)


class CarboxylicAcid(
    _generic_functional_group.GenericFunctionalGroup,
):
    """
    Represents a carboxylic acid functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][carbon](=[oxygen1])[oxygen2][hydrogen]``.

    """

    def __init__(
        self,
        carbon: _atoms.C,
        # O is not an ambiguous name.
        oxygen1: _atoms.O,  # noqa
        # O is not an ambiguous name.
        oxygen2: _atoms.O,  # noqa
        hydrogen: _atoms.H,
        atom: _atoms.Atom,
        bonders: tuple[_atoms.Atom, ...],
        deleters: tuple[_atoms.Atom, ...],
        placers: typing.Optional[tuple[_atoms.Atom, ...]] = None,
    ):
        """
        Initialize a :class:`.CarboxylicAcid` instance.

        Parameters:

            carbon:
                The ``[carbon]`` atom.

            oxygen1:
                The ``[oxygen1]`` atom.

            oxygen2:
                The ``[oxygen2]`` atom.

            hydrogen:
                The ``[hydrogen]`` atom.

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
            atoms=(carbon, oxygen1, oxygen2, hydrogen, atom),
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )
        self._carbon = carbon
        self._oxygen1 = oxygen1
        self._oxygen2 = oxygen2
        self._hydrogen = hydrogen
        self._atom = atom

    def get_carbon(self) -> _atoms.C:
        """
        Get the ``[carbon]`` atom.

        Returns:

            The ``[carbon]`` atom.

        """

        return self._carbon

    # O is not an ambiguous name.
    def get_oxygen1(self) -> _atoms.O:  # noqa
        """
        Get the ``[oxygen1]`` atom.

        Returns:

            The ``[oxygen]`` atom.

        """

        return self._oxygen1

    # O is not an ambiguous name.
    def get_oxygen2(self):  # noqa
        """
        Get the ``[oxygen2]`` atom.

        Returns:

            The ``[oxygen2]`` atom.

        """

        return self._oxygen2

    def get_hydrogen(self) -> _atoms.H:
        """
        Get the ``[hydrogen]`` atom.

        Returns:

            The ``[hydrogen]`` atom.

        """

        return self._hydrogen

    def get_atom(self) -> _atoms.Atom:
        """
        Get the ``[atom]`` atom.

        Returns:

            The ``[atom]`` atom.

        """

        return self._atom

    def clone(self) -> CarboxylicAcid:
        clone = self._clone()
        clone._carbon = self._carbon
        clone._oxygen1 = self._oxygen1
        clone._oxygen2 = self._oxygen2
        clone._hydrogen = self._hydrogen
        clone._atom = self._atom
        return clone

    def with_ids(
        self,
        id_map: dict[int, int],
    ) -> CarboxylicAcid:

        atom_map = _utilities.get_atom_map(
            id_map=id_map,
            atoms=(
                *self._atoms,
                *self._placers,
                *self._core_atoms,
                *self._bonders,
                *self._deleters,
                self._carbon,
                self._oxygen1,
                self._oxygen2,
                self._hydrogen,
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
        clone._oxygen1 = atom_map.get(
            self._oxygen1.get_id(),
            self._oxygen1,
        )
        clone._oxygen2 = atom_map.get(
            self._oxygen2.get_id(),
            self._oxygen2,
        )
        clone._hydrogen = atom_map.get(
            self._hydrogen.get_id(),
            self._hydrogen,
        )
        clone._atom = atom_map.get(
            self._atom.get_id(),
            self._atom,
        )
        return clone

    def __repr__(self) -> str:
        return (
            f'{self.__class__.__name__}('
            f'{self._carbon}, {self._oxygen1}, {self._oxygen2}, '
            f'{self._hydrogen}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
