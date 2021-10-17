"""
Atom Batch
==========

"""

import numpy as np
from collections import abc

from ....reactions import NewAtom
from ....atom import Atom
from ....atom_info import AtomInfo


__all__ = (
    'AtomBatch',
)


class AtomBatch:
    """
    A batch of atoms.

    """

    _atoms: list[Atom]
    _positions: list[np.ndarray]
    _atom_infos: list[AtomInfo]
    _atom_map: dict[int, Atom]

    __slots__ = ['_atoms', '_atom_infos', '_atom_map', '_positions']

    def __init__(
        self,
        atoms: abc.Iterable[NewAtom],
        num_atoms: int,
    ) -> None:
        """
        Initialize an :class:`.AtomBatch` instance.

        Parameters:

            atoms:
                The atoms, which should be added to the batch.

            num_atoms:
                The number of atoms in the molecule being constructed,
                before atoms in this batch are taken into account.

        """

        _atoms = self._atoms = []
        positions = self._positions = []
        atom_infos = self._atom_infos = []
        atom_map = self._atom_map = {}

        for id_, new_atom in enumerate(atoms, num_atoms):
            atom = new_atom.get_atom()
            position = new_atom.get_position()

            _atoms.append(atom.with_id(id_))
            atom_map[atom.get_id()] = _atoms[-1]
            atom_infos.append(AtomInfo(
                atom=_atoms[-1],
                building_block_atom=None,
                building_block=None,
                building_block_id=None,
            ))
            positions.append(position)

    def get_positions(self) -> abc.Iterator[np.ndarray]:
        """
        Yield the positions of atoms in the batch.

        Yields:

            The position of an atom in the batch.

        """

        yield from self._positions

    def get_atoms(self) -> abc.Iterator[Atom]:
        """
        Yield the atoms in the batch.

        Yields:

            An atom.

        """

        yield from self._atoms

    def get_atom_infos(self) -> abc.Iterator[AtomInfo]:
        """
        Yield info about the atoms in the batch.

        Yields:

            Info about an atom in the batch.

        """

        yield from self._atom_infos

    def get_atom_map(self) -> dict[int, Atom]:
        """
        Get a mapping from the old atom id to the new atom.

        Returns:

            Maps the id of an atom provided to the initializer, to
            the new atom held by the batch, which has an updated id.

        """

        return dict(self._atom_map)
