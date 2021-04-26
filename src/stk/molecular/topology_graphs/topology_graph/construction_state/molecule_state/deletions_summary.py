"""
Deletions Summary
=================

"""

import numpy as np

from .....atoms import AtomInfo
from .....bonds import BondInfo


class _DeletionsSummary:
    """
    A summary of deletion results.

    """

    __slots__ = [
        '_atoms',
        '_atom_infos',
        '_bonds',
        '_bond_infos',
        '_position_matrix',
        '_deleted_atom_ids',
        '_deleted_bond_ids',

        '_valid_atoms',
        '_valid_atom_infos',
        '_valid_bonds',
        '_valid_bond_infos',
        '_valid_positions',
    ]

    def __init__(
        self,
        atoms,
        atom_infos,
        bonds,
        bond_infos,
        position_matrix,
        deleted_atom_ids,
        deleted_bond_ids,
    ):
        """
        Initialize a :class:`._DeletionsSummary` instance.

        Parameters
        ----------
        atoms : :class:`iterable` of :class:`.Atom`
            Atoms, some of which are to be deleted.

        atom_infos : :class:`iterable` of :class:`.AtomInfo`
            Info on every atom in `atoms`.

        bonds : :class:`iterable` of :class:`.Bond`
            The bonds of the molecule being constructed.

        bond_infos : :class:`iterable` of :class:`.BondInfo`
            Info on every bond in `bonds`.

        position_matrix : :class:`numpy.ndarray`
            The position matrix for all the `atoms`.

        deleted_atom_ids : :class:`iterable` of :class:`int`
            The ids of `atoms`, which should be deleted.

        deleted_bond_ids : :class:`iterable` of :class:`._BondId`
            Ids of bonds that should be deleted.

        """

        self._atoms = tuple(atoms)
        self._atom_infos = tuple(atom_infos)
        self._bonds = tuple(bonds)
        self._bond_infos = tuple(bond_infos)
        self._position_matrix = np.array(position_matrix)
        self._deleted_atom_ids = set(deleted_atom_ids)
        self._deleted_bond_ids = set(deleted_bond_ids)

        self._valid_atoms = []
        self._valid_atom_infos = []
        self._valid_bonds = []
        self._valid_bond_infos = []
        self._valid_positions = []
        self._with_valid_data()

    def _with_valid_data(self):
        """
        Add the data of atoms which aren't deleted.

        Returns
        -------
        None : :class:`NoneType`

        """

        atoms = self._atoms
        atom_infos = self._atom_infos
        position_matrix = self._position_matrix
        deleted_atom_ids = self._deleted_atom_ids
        deleted_bond_ids = self._deleted_bond_ids

        def valid_atom(atom):
            return atom.get_id() not in deleted_atom_ids

        valid_atoms = self._valid_atoms
        valid_atom_infos = self._valid_atom_infos
        valid_positions = self._valid_positions
        atom_map = {}

        def with_atom(atom):
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

        def valid_bond(bond_data):
            index, bond = bond_data
            atom1_id = bond.get_atom1().get_id()
            atom2_id = bond.get_atom2().get_id()
            return (
                atom1_id not in deleted_atom_ids
                and atom2_id not in deleted_atom_ids
                and (atom1_id, atom2_id) not in deleted_bond_ids
            )

        bonds = self._bonds
        bond_infos = self._bond_infos
        valid_bonds = self._valid_bonds
        valid_bond_infos = self._valid_bond_infos

        def with_bond(index, bond):
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

    def get_atoms(self):
        """
        Yield the atoms in the summary.

        Yields
        ------
        :class:`.Atom`
            An atom.

        """

        yield from self._valid_atoms

    def get_atom_infos(self):
        """
        Yield infos about atoms in the summary.

        Yields
        ------
        :class:`.AtomInfo`
            Info about an atom.

        """

        yield from self._valid_atom_infos

    def get_bonds(self):
        """
        Yield the bonds in the summary.

        Yields
        ------
        :class:`.Bond`
            A bond.

        """

        yield from self._valid_bonds

    def get_bond_infos(self):
        """
        Yield infos about the bonds in the summary.

        Yields
        ------
        :class:`.BondInfo`
            Info about a bond.

        """

        yield from self._valid_bond_infos

    def get_positions(self):
        """
        Yield the positions of atoms held by the summary.

        Yields
        ------
        :class:`numpy.ndarray`
            The position of an atom held by the summary.

        """

        yield from self._valid_positions
