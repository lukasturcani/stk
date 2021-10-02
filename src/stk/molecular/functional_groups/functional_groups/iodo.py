"""
Iodo
====

"""

from __future__ import annotations

import typing

from . import utilities as _utilities
from . import generic_functional_group as _generic_functional_group
from ... import atoms as _atoms


__all__ = (
    'Iodo',
)


class Iodo(
    _generic_functional_group.GenericFunctionalGroup,
):
    """
    Represents an iodo functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[iodine][atom]``.

    """

    def __init__(
        self,
        # I is not an ambiguous name.
        iodine: _atoms.I,  # noqa
        atom: _atoms.Atom,
        bonders: tuple[_atoms.Atom, ...],
        deleters: tuple[_atoms.Atom, ...],
        placers: typing.Optional[tuple[_atoms.Atom, ...]] = None,
    ) -> None:
        """
        Initialize a :class:`.Iodo` instance.

        Parameters:

            iodine:
                The ``[iodine]`` atom.

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
            atoms=(iodine, atom),
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )
        self._iodine = iodine
        self._atom = atom

    # I is not an ambiguous name.
    def get_iodine(self) -> _atoms.I:  # noqa
        """
        Get the ``[iodine]`` atom.

        Returns:

            The ``[iodine]`` atom.

        """

        return self._iodine

    def get_atom(self) -> _atoms.Atom:
        """
        Get the ``[atom]`` atom.

        Returns:

            The ``[atom]`` atom.

        """

        return self._atom

    def clone(self) -> Iodo:
        clone = self._clone()
        clone._iodine = self._iodine
        clone._atom = self._atom
        return clone

    def with_ids(
        self,
        id_map: dict[int, int],
    ) -> Iodo:

        atom_map = _utilities.get_atom_map(
            id_map=id_map,
            atoms=(
                *self._atoms,
                *self._placers,
                *self._core_atoms,
                *self._bonders,
                *self._deleters,
                self._iodine,
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
        clone._iodine = atom_map.get(
            self._iodine.get_id(),
            self._iodine,
        )
        clone._atom = atom_map.get(
            self._atom.get_id(),
            self._atom,
        )
        return clone

    def __repr__(self) -> str:
        return (
            f'{self.__class__.__name__}('
            f'{self._iodine}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
