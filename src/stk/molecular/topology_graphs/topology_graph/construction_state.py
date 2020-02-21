import numpy as np
from collections import defaultdict

from ...atom_info import AtomInfo


class ConstructionState:
    """

    """

    def __init__(self, building_block_vertices, edges, vertex_edges):
        """

        """

        self._vertex_building_blocks = {
            vertex.get_id(): building_block
            for building_block, vertices
            in building_block_vertices.items()
            for vertex in vertices
        }
        self._vertices = {
                vertex.get_id(): vertex
                for vertices in building_block_vertices.values()
                for vertex in vertices
        }
        self._vertex_edges = dict(vertex_edges)
        self._edges = edges
        self._position_matrx = []
        self._atoms = []
        self._atom_infos = []
        self._bonds = []
        self._building_block_counts = defaultdict(int)
        self._reaction_infos = []
        self._edge_functional_groups = defaultdict(list)

    def clone(self):
        clone = self.__class__.__new__(self.__class__)
        clone._vertex_building_blocks = dict(
            self._vertex_building_blocks
        )
        clone._vertices = dict(self._vertices)
        clone._vertex_edges = dict(self._vertex_edges)
        clone._edges = self._edges
        clone._position_matrix = list(self._position_matrix)
        clone._atoms = list(self._atoms)
        clone._atom_infos = list(self._atom_infos)
        clone._bonds = list(self._bonds)
        clone._building_block_counts = defaultdict.copy(
            self._building_block_counts
        )
        clone._reaction_infos = list(self._reaction_infos)
        clone._edge_functional_groups = defaultdict.copy(
            self._edge_functional_groups
        )

    def _with_placement_results(self, building_blocks, results):
        """
        Modify the state.

        """

        results_ = zip(building_blocks, results)
        for index, (building_block, result) in enumerate(results_):
            self._position_matrix.extend(result.position_matrix)
            self._building_block_counts[building_block] += 1

            # atom_map maps atoms in the original building block to
            # the new atoms the constructed molecule should hold.
            # This means their atom_ids are updated.
            atom_map = {}

            # Create atom_map, add the new atoms and create AtomInfos.
            for atom in building_block.get_atoms():
                new_atom = atom.with_id(len(self._atoms))
                atom_map[atom.get_id()] = new_atom

                self._atoms.append(atom)
                self._atom_infos.append(
                    AtomInfo(
                        atom=new_atom,
                        building_block=building_block,
                        building_block_index=index,
                    )
                )

            # Add bonds holding the new atoms.
            for bond in building_block.get_bonds():
                self._bonds.append(bond.with_atoms(atom_map))

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

        return self._vertex[vertex_id]

    def get_edge(self, edge_id):
        return self._edges[edge_id]

    def get_edges(self, vertex_id):
        return self._vertex_edges[vertex_id]

    def get_edge_group_functional_groups(self, edge_group):
        for edge_id in edge_group.get_edge_ids():
            yield from self._edge_functional_groups[edge_id]

    def with_reaction_results(self, results):
        pass

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

    def get_reaction_infos(self):
        yield from self._reaction_infos

    def get_building_block_counts(self):
        return dict(self._building_block_counts)
