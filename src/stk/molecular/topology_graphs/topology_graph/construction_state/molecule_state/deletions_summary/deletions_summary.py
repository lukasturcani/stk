import numpy as np

from .....atoms import AtomInfo
from .....bonds import BondInfo


class _DeletionsSummary:
    def __init__(
        self,
        atoms,
        atom_infos,
        bonds,
        bond_infos,
        position_matrix,
        deleted_ids,
    ):
        self._atoms = tuple(atoms)
        self._atom_infos = tuple(atom_infos)
        self._bonds = tuple(bonds)
        self._bond_infos = tuple(bond_infos)
        self._position_matrix = np.array(position_matrix)
        self._deleted_ids = set(deleted_ids)

        self._valid_atoms = []
        self._valid_atom_infos = []
        self._valid_bonds = []
        self._valid_bond_infos = []
        self._valid_positions = []
        self._with_valid_data()

    def _with_valid_data(self):
        atoms = self._atoms
        atom_infos = self._atom_infos
        position_matrix = self._position_matrix
        deleted_ids = self._deleted_ids

        def valid_atom(atom):
            return atom.get_id() not in deleted_ids

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
                    building_block=info.get_building_block(),
                    building_block_id=info.get_building_block_id(),
                )
            )
            valid_positions.append(position_matrix[atom_id])

        for atom in filter(valid_atom, atoms):
            with_atom(atom)

        def valid_bond(bond_data):
            index, bond = bond_data
            return (
                bond.get_atom1().get_id() not in deleted_ids
                and bond.get_atom2().get_id() not in deleted_ids
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
        yield from self._valid_atoms

    def get_atom_infos(self):
        yield from self._valid_atom_ids

    def get_bonds(self):
        yield from self._valid_bonds

    def get_bond_infos(self):
        yield from self._valid_bond_infos

    def get_positions(self):
        yield from self._valid_positions
