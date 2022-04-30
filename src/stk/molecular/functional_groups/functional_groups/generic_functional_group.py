"""
Generic Functional Group
========================

"""


from __future__ import annotations

import typing

from ...atoms import Atom
from .functional_group import FunctionalGroup

_T = typing.TypeVar('_T', bound='GenericFunctionalGroup')


class GenericFunctionalGroup(FunctionalGroup):
    """
    A functional group which defines general atomic classes.

    *Bonders* are atoms which should have bonds added by a
    :class:`.Reaction`. *Deleters* are atoms which should be removed
    by a :class:`.Reaction`.

    This interface allows the same reactions to be carried out across
    different functional groups, without knowing which specific
    elements it holds.

    """

    def __init__(
        self,
        atoms: tuple[Atom, ...],
        bonders: tuple[Atom, ...],
        deleters: tuple[Atom, ...],
        placers: typing.Optional[tuple[Atom, ...]] = None,
    ) -> None:
        """
        Initialize a :class:`.GenericFunctionalGroup`.

        Parameters:

            atoms:
                The atoms in the functional group.

            bonders:
                The bonder atoms in the functional group.

            deleters:
                The deleter atoms in the functional group.

            placers:
                The placer atoms of the functional group. If ``None``,
                the `bonders` will be used.

        """

        deleter_set = set(atom.get_id() for atom in deleters)
        super().__init__(
            atoms=atoms,
            placers=bonders if placers is None else placers,
            core_atoms=tuple(
                atom for atom in atoms
                if atom.get_id() not in deleter_set
            ),
        )
        self._bonders = bonders
        self._deleters = deleters

    def _clone(self: _T) -> _T:
        clone = super()._clone()
        clone._bonders = self._bonders
        clone._deleters = self._deleters
        return clone

    def clone(self) -> GenericFunctionalGroup:
        return self._clone()

    def with_atoms(
        self,
        atom_map: dict[int, Atom],
    ) -> GenericFunctionalGroup:

        return GenericFunctionalGroup(
            atoms=tuple(
                atom_map.get(a.get_id(), a) for a in self._atoms
            ),
            bonders=tuple(
                atom_map.get(a.get_id(), a) for a in self._bonders
            ),
            deleters=tuple(
                atom_map.get(a.get_id(), a) for a in self._deleters
            ),
            placers=tuple(
                atom_map.get(a.get_id(), a) for a in self._placers
            ),
        )

    def _with_ids(self: _T, id_map: dict[int, int]) -> _T:
        super()._with_ids(id_map)
        self._bonders = tuple(
            bonder.with_id(
                id=id_map.get(
                    bonder.get_id(),
                    bonder.get_id(),
                )
            ) for bonder in self._bonders
        )
        self._deleters = tuple(
            deleter.with_id(
                id=id_map.get(
                    deleter.get_id(),
                    deleter.get_id(),
                ),
            ) for deleter in self._deleters
        )
        return self

    def with_ids(
        self,
        id_map: dict[int, int],
    ) -> GenericFunctionalGroup:

        return self.clone()._with_ids(id_map)

    def get_bonders(self) -> typing.Iterator[Atom]:
        """
        Yield bonder atoms in the functional group.

        These are atoms which have bonds added during
        :class:`.ConstructedMolecule` construction.

        Yields:
            A bonder atom.

        """

        yield from self._bonders

    def get_num_bonders(self) -> int:
        """
        Get the number of bonder atoms.

        Returns:
            The number of bonder atoms.

        """

        return len(self._bonders)

    def get_bonder_ids(self) -> typing.Iterator[int]:
        """
        Yield the ids of bonder atoms.

        Yields:
            The id of a bonder :class:`.Atom`.

        """

        yield from (a.get_id() for a in self._bonders)

    def get_deleters(self) -> typing.Iterator[Atom]:
        """
        Yield the deleter atoms in the functional group.

        These are atoms which are removed during
        :class:`.ConstructedMolecule` construction.

        Yields:
            A deleter atom.

        """

        yield from self._deleters

    def get_deleter_ids(self) -> typing.Iterator[int]:
        """
        Yield the ids of deleter atoms.

        Yields:
            The id of a deleter :class:`.Atom`.

        """

        yield from (a.get_id() for a in self._deleters)

    def __repr__(self) -> str:
        return (
            f'{self.__class__.__name__}('
            f'atoms={self._atoms}, '
            f'bonders={self._bonders}, '
            f'deleters={self._deleters}'
            ')'
        )
