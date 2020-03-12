import numpy as np

from .....atoms import AtomInfo
from .....bonds import BondInfo


class _PlacementsSummary:
    def __init__(self, building_blocks, placement_results, next_id):
        self._atoms = []
        self._atom_infos = []
        self._bonds = []
        self._bond_infos = []
        self._edge_functional_groups = {}
        self._position_matrices = []

        results = zip(building_blocks, placement_results)
        for id_, (building_block, result) in enumerate(results, next_id):
            self._with_placement_result(building_block, id_, result)

    def _with_placement_result(
        self,
        building_block,
        building_block_id,
        result,
    ):
        atoms = self._atoms
        atom_infos = self._atom_infos
        atom_map = {}

        def with_atom(atom):
            atoms.append(atom.with_id(len(atoms)))
            atom_map[atom.get_id()] = atoms[-1]
            atom_infos.append(
                AtomInfo(
                    atom=atoms[-1],
                    building_block=building_block,
                    building_block_id=building_block_id,
                )
            )

        for atom in building_block.get_atoms():
            with_atom(atom)

        bonds = self._bonds
        bond_infos = self._bond_infos

        def with_bond(bond):
            bonds.append(bond.with_atoms(atom_map))
            bond_infos.append(
                BondInfo(
                    bond=bonds[-1],
                    building_block=building_block,
                    building_block_id=building_block_id,
                )
            )
        for bond in building_block.get_bonds():
            with_bond(bond)

        self._position_matrices.append(result.position_matrix)

    def get_atoms(self):
        yield from self._atoms

    def get_atom_infos(self):
        yield from self._atom_infos

    def get_bonds(self):
        yield from self._bonds

    def get_bond_infos(self):
        yield from self._bond_infos

    def get_position_matrix(self):
        return np.vstack(self._position_matrices)

    def get_edge_functional_groups(self):
        yield from self._edge_functional_groups.items()

    def _with_placement_results(
        self,
        vertices,
        edges,
        building_blocks,
        results,
    ):
        """
        Modify the state.

        """

        # Doing a vstack after the loop should be faster than doing one
        # within the loop.
        position_matrices = [self._position_matrix]
        results_ = zip(building_blocks, results)
        for index, (building_block, result) in enumerate(results_):
            position_matrices.append(result.position_matrix)
            atom_map = self._add_atoms(index, building_block)
            self._add_bonds(index, building_block, atom_map)
            self._add_edge_functional_groups(
                building_block=building_block,
                result=result,
                atom_map=atom_map,
            )
        self._position_matrix = np.vstack(position_matrices)
        return self

    def _add_atoms(self, index, building_block):
        atom_map = {}
        for atom in building_block.get_atoms():
            new_atom = atom.with_id(len(self._atoms))
            atom_map[atom.get_id()] = new_atom
            self._atoms.append(new_atom)
            self._atom_infos.append(
                AtomInfo(
                    atom=new_atom,
                    building_block=building_block,
                    building_block_id=index,
                )
            )
        return atom_map

    def _add_bonds(self, index, building_block, atom_map):
        for bond in building_block.get_bonds():
            new_bond = bond.with_atoms(atom_map)
            self._bonds.append(new_bond)
            self._bond_infos.append(
                BondInfo(
                    bond=new_bond,
                    building_block=building_block,
                    building_block_id=index,
                )
            )

    def _add_edge_functional_groups(
        self,
        building_block,
        result,
        atom_map,
    ):
        # Add edge to functional group mappings.
        functional_groups = building_block.get_functional_groups(
            fg_ids=result.functional_group_edges,
        )
        edge_ids = result.functional_group_edges.values()
        functional_group_edges = zip(functional_groups, edge_ids)
        for functional_group, edge_id in functional_group_edges:
            self._edge_functional_groups[edge_id].append(
                functional_group.with_atoms(atom_map)
            )
