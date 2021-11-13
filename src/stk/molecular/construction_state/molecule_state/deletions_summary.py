"""
Deletions Summary
=================

"""

from collections import abc
import numpy as np

from .reactions_summary import BondId
from ...atom import Atom
from ...atom_info import AtomInfo
from ...bond import Bond
from ...bond_info import BondInfo


__all__ = (
    'DeletionsSummary',
)


class DeletionsSummary:
    """
    A summary of deletion results.

    """

    _valid_atoms: list[Atom]
    _valid_atom_infos: list[AtomInfo]
    _valid_bonds: list[Bond]
    _valid_bond_infos: list[BondInfo]
    _valid_positions: list[np.ndarray]

    __slots__ = [
        '_valid_atoms',
        '_valid_atom_infos',
        '_valid_bonds',
        '_valid_bond_infos',
        '_valid_positions',
    ]

    def __init__(
        self,
        atoms: abc.Iterable[Atom],
        atom_infos: abc.Sequence[AtomInfo],
        bonds: abc.Iterable[Bond],
        bond_infos: abc.Sequence[BondInfo],
        position_matrix: np.ndarray,
        deleted_atom_ids: frozenset[int],
        deleted_bond_ids: frozenset[BondId],
    ) -> None:
        """
        Initialize a :class:`.DeletionsSummary` instance.

        Parameters:

            atoms:
                Atoms, some of which are to be deleted.

            atom_infos:
                Info on every atom in `atoms`.

            bonds:
                The bonds of the molecule being constructed.

            bond_infos:
                Info on every bond in `bonds`.

            position_matrix:
                The position matrix for all the `atoms`.

            deleted_atom_ids:
                The ids of `atoms`, which should be deleted.

            deleted_bond_ids:
                Ids of bonds that should be deleted.

        """

        self._valid_atoms = []
        self._valid_atom_infos = []
        self._valid_bonds = []
        self._valid_bond_infos = []
        self._valid_positions = []

        def valid_atom(atom: Atom) -> bool:
            return atom.get_id() not in deleted_atom_ids

        valid_atoms = self._valid_atoms
        valid_atom_infos = self._valid_atom_infos
        valid_positions = self._valid_positions
        atom_map = {}

        def with_atom(atom: Atom) -> None:
            atom_id = atom.get_id()
            valid_atoms.append(atom.with_id(len(valid_atoms)))
            atom_map[atom_id] = valid_atoms[-1]
            info = atom_infos[atom_id]
            valid_atom_infos.append(
                AtomInfo(
                    atom=valid_atoms[-1],
                    building_block_atom=info.get_building_block_atom(),
                    building_block=info.get_building_block(),
                    building_block_id=info.get_building_block_id(),
                )
            )
            valid_positions.append(position_matrix[atom_id])

        for atom in filter(valid_atom, atoms):
            with_atom(atom)

        def valid_bond(bond_data: tuple[int, Bond]) -> bool:
            index, bond = bond_data
            atom1_id = bond.get_atom1().get_id()
            atom2_id = bond.get_atom2().get_id()
            return (
                atom1_id not in deleted_atom_ids
                and atom2_id not in deleted_atom_ids
                and BondId(atom1_id, atom2_id) not in deleted_bond_ids
            )

        valid_bonds = self._valid_bonds
        valid_bond_infos = self._valid_bond_infos

        def with_bond(index: int, bond: Bond) -> None:
            valid_bonds.append(bond.with_atoms(atom_map))
            info = bond_infos[index]
            valid_bond_infos.append(
                BondInfo(
                    bond=valid_bonds[-1],
                    building_block=info.get_building_block(),
                    building_block_id=info.get_building_block_id(),
                )
            )

        for index, bond in filter(valid_bond, enumerate(bonds)):
            with_bond(index, bond)

    def get_atoms(self) -> abc.Iterator[Atom]:
        """
        Yield the atoms in the summary.

        Yields:

            An atom.

        """

        yield from self._valid_atoms

    def get_atom_infos(self) -> abc.Iterator[AtomInfo]:
        """
        Yield infos about atoms in the summary.

        Yields:

            Info about an atom.

        """

        yield from self._valid_atom_infos

    def get_bonds(self) -> abc.Iterator[Bond]:
        """
        Yield the bonds in the summary.

        Yields:

            A bond.

        """

        yield from self._valid_bonds

    def get_bond_infos(self) -> abc.Iterator[BondInfo]:
        """
        Yield infos about the bonds in the summary.

        Yields:

            Info about a bond.

        """

        yield from self._valid_bond_infos

    def get_positions(self) -> abc.Iterator[np.ndarray]:
        """
        Yield the positions of atoms held by the summary.

        Yields:

            The position of an atom held by the summary.

        """

        yield from self._valid_positions
