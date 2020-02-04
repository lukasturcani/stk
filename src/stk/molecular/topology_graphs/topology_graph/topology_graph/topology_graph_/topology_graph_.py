import numpy as np
import pathos
from collections import namedtuple

from .implementations import _Parallel, _Serial
from ...utilities import vector_angle


class TopologyGraph_:
    """
    A partial implementation of :class:`.TopologyGraph`.

    """

    def __init__(
        self,
        vertex_data,
        edge_data,
        construction_stages,
        num_processes,
    ):
        """
        Initialize an instance of :class:`.TopologyGraph`.

        Parameters
        ----------
        vertex_data : :class:`tuple` of :class:`.VertexData`
            The vertices which make up the graph.

        edge_data : :class:`tuple` of :class:`.EdgeData`
            The edges which make up the graph.

        construction_stages : :class:`tuple` of :class:`callable`
            A collection of callables, each of which takes a
            :class:`.Vertex` and returns ``True`` or ``False``.
            If the first :class:`callable` is applied to a  vertex in
            `vertices`, that vertex is is part of the first
            construction stage. The second :class:`callable` is then
            applied to all vertices not in the first stage and those
            which return ``True`` belong to the second stage and
            so on.

            Vertices which belong to the same construction stage
            all place building blocks together in parallel, before
            placement is done by any vertices which are part of a later
            stage. This breaks down parallel construction into
            serial stages if synchronization between stages is needed.

            If the topology graph is performing construction serially,
            then all vertices which belong to an earlier stage will
            place their building block before those at a later stage.

        num_processes : :class:`int`
            The number of parallel processes to create during
            :meth:`construct`.

        """

        vertex_data = tuple(self._with_ids(vertex_data))
        edge_data = tuple(self._with_ids(edge_data))
        self._vertices = tuple(
            data.get_vertex() for data in vertex_data
        )
        self._edges = tuple(
            data.get_edge() for data in edge_data
        )
        self._construction_stages = construction_stages

        if num_processes == 1:
            self._implementation = _Serial(
                stages=self._get_stages(),
                after_placement_stage=self._after_placement_stage,
            )
        else:
            self._implementation = _Parallel(
                stages=self._get_stages(),
                num_processes=num_processes,
                after_placement_stage=self._after_placement_stage,
            )

    @staticmethod
    def _with_ids(objects):
        for id, object in enumerate(objects):
            yield object.with_id(id)

    def _get_stages(self):
        stages = tuple(
            [] for i in range(len(self._construction_stages)+1)
        )
        for vertex in self._vertices:
            placed = False
            for i, stage in enumerate(self._construction_stages):
                if stage(vertex):
                    self._stages[i].append(vertex.get_id())
                    placed = True
                    break
            if not placed:
                self._stages[-1].append(vertex.get_id())
        return stages

    def _get_vertex_clones(self, mol, scale):
        for vertex in self.vertices:
            yield (
                vertex
                .clone()
                .set_constructed_molecule(mol)
                .apply_scale(scale)
            )

    def _get_edge_clones(self, scale):
        for edge in self.edges:
            yield edge.clone().apply_scale(scale)

    def _before_react(self, mol, vertices, edges):
        return vertices, edges

    def _prepare(self, mol):
        return

    def _place_building_blocks(self, state):
        return self._implementation._place_building_blocks(state)

    def _get_atom_map(self, mol, bb, bb_id):
        atom_map = {}
        for atom in bb.atoms:
            atom_clone = atom.clone()
            atom_clone.id = len(mol.atoms)
            atom_clone.building_block = bb
            atom_clone.building_block_id = bb_id
            atom_map[atom] = atom_clone
            mol.atoms.append(atom_clone)
        return atom_map

    def _assign_func_groups_to_edges(
        self,
        mol,
        bb,
        bb_id,
        edges,
        assignments
    ):
        atom_map = self._get_atom_map(mol, bb, bb_id)
        mol.func_groups.extend(
            fg.clone(atom_map) for fg in bb.func_groups
        )
        num_fgs = len(bb.func_groups)
        for fg_id, edge_id in assignments.items():
            edges[edge_id].assign_func_group(
                func_group=mol.func_groups[-num_fgs+fg_id]
            )
        return atom_map

    def _clean_up(self, mol):
        mol._position_matrix = np.array(mol._position_matrix).T
        for i, atom in enumerate(mol.atoms):
            atom.id = i

    def get_adjacency_list(self):
        raise NotImplementedError()

    def __str__(self):
        return repr(self)

    def __repr__(self):
        raise NotImplementedError()
