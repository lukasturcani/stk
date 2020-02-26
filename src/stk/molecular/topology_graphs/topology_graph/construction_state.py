import numpy as np
from collections import defaultdict

from ...atom_info import AtomInfo
from ...bond_info import BondInfo


class ConstructionState:
    """

    """

    def __init__(
        self,
        building_block_vertices,
        edges,
        vertex_edges,
        scale,
    ):
        """

        """

        self._vertex_building_blocks = {
            vertex.get_id(): building_block
            for building_block, vertices
            in building_block_vertices.items()
            for vertex in vertices
        }
        self._vertices = {
                vertex.get_id(): vertex.with_scale(scale)
                for vertices in building_block_vertices.values()
                for vertex in vertices
        }
        self._vertex_edges = {
            vertex_id: tuple(
                edge.with_scale(scale) for edge in edges
            )
            for vertex_id, edges in vertex_edges.items()
        }
        self._edges = tuple(edge.with_scale(scale) for edge in edges)
        self._position_matrix = np.empty((0, 3), dtype=np.float64)
        self._atoms = []
        self._atom_infos = []
        self._bonds = []
        self._bond_infos = []
        self._building_block_counts = defaultdict(int)
        self._edge_functional_groups = defaultdict(list)

    def clone(self):
        clone = self.__class__.__new__(self.__class__)
        clone._vertex_building_blocks = dict(
            self._vertex_building_blocks
        )
        clone._vertices = dict(self._vertices)
        clone._vertex_edges = dict(self._vertex_edges)
        clone._edges = self._edges
        clone._position_matrix = np.array(self._position_matrix)
        clone._atoms = list(self._atoms)
        clone._atom_infos = list(self._atom_infos)
        clone._bonds = list(self._bonds)
        clone._bond_infos = list(self._bond_infos)
        clone._building_block_counts = defaultdict.copy(
            self._building_block_counts
        )
        clone._edge_functional_groups = defaultdict.copy(
            self._edge_functional_groups
        )
        return clone

    def _with_placement_results(self, building_blocks, results):
        """
        Modify the state.

        """

        # Doing a vstack after the loop should be faster than doing one
        # within the loop.
        position_matrices = [self._position_matrix]
        results_ = zip(building_blocks, results)
        for index, (building_block, result) in enumerate(results_):
            position_matrices.append(result.position_matrix)
            self._building_block_counts[building_block] += 1

            # atom_map maps atoms in the original building block to
            # the new atoms the constructed molecule should hold.
            # This means their atom_ids are updated.
            atom_map = {}

            # Create atom_map, add the new atoms and create AtomInfos.
            for atom in building_block.get_atoms():
                new_atom = atom.with_id(len(self._atoms))
                atom_map[atom.get_id()] = new_atom

                self._atoms.append(new_atom)
                self._atom_infos.append(
                    AtomInfo(
                        atom=new_atom,
                        building_block=building_block,
                        building_block_index=index,
                    )
                )

            # Add bonds holding the new atoms.
            for bond in building_block.get_bonds():
                new_bond = bond.with_atoms(atom_map)
                self._bonds.append(new_bond)
                self._bond_infos.append(
                    BondInfo(
                        bond=new_bond,
                        building_block=building_block,
                        building_block_index=index,
                    )
                )

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

        self._position_matrix = np.vstack(position_matrices)
        return self

    def with_placement_results(self, building_blocks, results):
        return self.clone()._with_placement_results(
            building_blocks=building_blocks,
            results=results,
        )

    def get_building_block(self, vertex_id):
        """

        """

        return self._vertex_building_blocks[vertex_id]

    def get_vertex(self, vertex_id):
        """

        """

        return self._vertices[vertex_id]

    def get_num_vertices(self):
        return len(self._vertices)

    def get_edge(self, edge_id):
        return self._edges[edge_id]

    def get_num_edges(self):
        return len(self._edges)

    def get_edges(self, vertex_id):
        return self._vertex_edges[vertex_id]

    def get_edge_group_functional_groups(self, edge_group):
        for edge_id in edge_group.get_edge_ids():
            yield from self._edge_functional_groups[edge_id]

    def _add_new_atoms_and_bonds(self, reactions, results):
        new_atom_positions = []
        deleted_atoms = set()
        for reaction, result in zip(reactions, results):
            deleted_atoms.update(
                atom.get_id() for atom in result.deleted_atoms
            )

            atom_map = {}
            for atom, position in result.new_atoms:
                new_atom = atom.with_id(len(self._atoms))
                atom_map[atom.get_id()] = new_atom
                self._atoms.append(new_atom)
                self._atom_infos.append(
                    AtomInfo(
                        atom=new_atom,
                        building_block=None,
                        building_block_index=None,
                    )
                )
                new_atom_positions.append(position)

            for bond in result.new_bonds:
                new_bond = bond.with_atoms(atom_map)
                self._bonds.append(new_bond)
                self._bond_infos.append(
                    BondInfo(
                        bond=new_bond,
                        building_block=None,
                        building_block_index=None,
                    )
                )
        if new_atom_positions:
            self._position_matrix = np.vstack([
                self._position_matrix,
                new_atom_positions,
            ])
        return deleted_atoms

    def _remove_deleted_atoms(self, deleted_atoms):
        atom_map = {}
        valid_atoms = filter(
            lambda atom: atom.get_id() not in deleted_atoms,
            self._atoms,
        )
        self._atoms = []
        position_matrix = []
        for atom in valid_atoms:
            new_atom = atom.with_id(len(self._atoms))
            atom_map[atom.get_id()] = new_atom
            self._atoms.append(new_atom)
            position_matrix.append(
                self._position_matrix[atom.get_id()]
            )
        self._position_matrix = np.array(position_matrix)

        valid_atom_infos = filter(
            lambda atom_info: atom_info.atom.get_id() in atom_map,
            self._atom_infos
        )
        self._atom_infos = [
            AtomInfo(
                atom=atom_map[atom_info.atom.get_id()],
                building_block=atom_info.building_block,
                building_block_index=atom_info.building_block_index,
            )
            for atom_info in valid_atom_infos
        ]
        return atom_map

    def _update_bonds(self, deleted_atoms, atom_map):
        valid_bonds = filter(
            lambda bond: (
                bond.get_atom1().get_id() not in deleted_atoms
                and bond.get_atom2().get_id() not in deleted_atoms
            ),
            self._bonds,
        )
        self._bonds = [
            bond.with_atoms(atom_map) for bond in valid_bonds
        ]

        valid_bond_infos = filter(
            lambda bond_info: (
                bond_info.bond.get_atom1().get_id()
                not in deleted_atoms
                and bond_info.bond.get_atom2().get_id()
                not in deleted_atoms
            ),
            self._bond_infos,
        )
        self._bond_infos = [
            BondInfo(
                bond=self._bonds[index],
                building_block=bond_info.building_block,
                building_block_index=bond_info.building_block_index,
            )
            for index, bond_info in enumerate(valid_bond_infos)
        ]

    def _with_reaction_results(self, reactions, results):
        deleted_atoms = self._add_new_atoms_and_bonds(
            reactions=reactions,
            results=results,
        )
        atom_map = self._remove_deleted_atoms(deleted_atoms)
        self._update_bonds(deleted_atoms, atom_map)
        return self

    def with_reaction_results(self, reactions, results):
        return self.clone()._with_reaction_results(reactions, results)

    def _with_vertices(self, vertices):
        self._vertices = {
            vertex.get_id(): vertex for vertex in vertices
        }
        return self

    def with_vertices(self, vertices):
        return self.clone()._with_vertices(vertices)

    def get_position_matrix(self):
        """

        """

        return np.array(self._position_matrix)

    def get_atoms(self):
        yield from self._atoms

    def get_bonds(self):
        yield from self._bonds

    def get_atom_infos(self):
        yield from self._atom_infos

    def get_bond_infos(self):
        yield from self._bond_infos

    def get_building_block_counts(self):
        return dict(self._building_block_counts)
