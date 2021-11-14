"""
Atom Batch
==========

"""

from __future__ import annotations

from typing import Iterable

from ......atoms import Atom, AtomInfo
from ......molecules import BuildingBlock


class _AtomBatch:
    """
    A batch of atoms.

    """

    __slots__ = ['_atoms', '_atom_infos', '_id_map']

    _atoms: tuple[Atom, ...]
    _atom_infos: tuple[AtomInfo, ...]
    _id_map: dict[int, int]

    def __init__(
        self,
        atoms: Iterable[Atom],
        num_atoms: int,
        building_block: BuildingBlock,
        building_block_id: int,
    ) -> None:
        """
        Initialize an :class:`._AtomBatch` instance.

        Parameters:

            atoms:
                The atoms, which should be added to the batch.

            num_atoms:
                The number of atoms in the molecule being constructed,
                before atoms in this batch are taken into account.

            building_block:
                The building block from which the atoms originate.

            building_block_id:
                An id, unique to that building block and placement.

        """

        _atoms = []
        atom_infos = []
        self._id_map = {}

        for id_, atom in enumerate(atoms, num_atoms):
            _atoms.append(atom.with_id(id_))
            self._id_map[atom.get_id()] = _atoms[-1].get_id()
            atom_infos.append(
                AtomInfo(
                    atom=_atoms[-1],
                    building_block_atom=atom,
                    building_block=building_block,
                    building_block_id=building_block_id,
                )
            )

        self._atoms = tuple(_atoms)
        self._atom_infos = tuple(atom_infos)

    def get_atoms(self) -> Iterable[Atom]:
        """
        Yield the atoms in the batch.

        Yields:

            An atom.

        """

        yield from self._atoms

    def get_atom_infos(self) -> Iterable[AtomInfo]:
        """
        Yield info about the atoms in the batch.

        Yields:

            Info about an atom in the batch.

        """

        yield from self._atom_infos

    def get_id_map(self) -> dict[int, int]:
        """
        Get a mapping from the old atom id to the new atom.

        Returns:

            Maps the id of an atom provided to the initializer, to
            the id of the new, corresponding, atom held by the batch.

        """

        return dict(self._id_map)
