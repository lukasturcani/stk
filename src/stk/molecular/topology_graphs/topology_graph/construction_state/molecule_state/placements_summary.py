import numpy as np
from collections import defaultdict

from .....atoms import AtomInfo
from .....bonds import BondInfo


class _PlacementsSummary:
    def __init__(self, building_blocks, placement_results, next_id):
        self._atoms = []
        self._atom_infos = []
        self._bonds = []
        self._bond_infos = []
        self._edge_functional_groups = defaultdict(list)
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

        def with_edge_functional_groups():
            functional_groups = building_block.get_functional_groups(
                fg_ids=result.functional_group_edges,
            )
            edge_ids = result.functional_group_edges.values()
            functional_group_edges = zip(functional_groups, edge_ids)
            for functional_group, edge_id in functional_group_edges:
                self._edge_functional_groups[edge_id].append(
                    functional_group.with_atoms(atom_map)
                )

        with_edge_functional_groups()

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
