import numpy as np
import pathos
from collections import namedtuple

from .implementations import Parallel, Serial
from ...utilities import vector_angle


class TopologyGraph:
    """
    Represents topology graphs of :class:`.ConstructedMolecule`.

    The topology graph is an abstract representation of a constructed
    molecule. The vertices indicate where building blocks are placed
    and the edges indicate which building blocks have bonds formed
    between them by the construction process.

    Vertices are responsible for placing the building block molecules.
    By initializing the vertices with different settings, they can
    position the building block molecules differently and therefore
    allow the user to easily specify a different structural isomer.

    Once a building block is placed on a vertex, the functional groups
    on the building block must be assigned to the different edges
    connected to the vertex. The number of functional groups in the
    building block must match the number of edges connected to the
    vertex.

    Once the functional groups are assigned to edges, each edge
    represents a reaction between the functional groups assigned to it.
    Note that an edge can be assigned more than two functional groups,
    in case you are dealing with something really exotic. The
    functional groups are then matched to an appropriate reaction,
    which generally creates bonds between the atoms of the functional
    groups. After this you will end up with a
    :class:`.ConstructedMolecule`.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    def __init__(
        self,
        vertex_data,
        edge_data,
        construction_stages,
        num_processes
    ):
        """
        Initialize an instance of :class:`.TopologyGraph`.

        Parameters
        ----------
        vertices : :class:`tuple` of :class:`.VertexData`
            The vertices which make up the graph.

        edges : :class:`tuple` of :class:`.EdgeData`
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

        self._set_data_ids(vertex_data)
        self._set_data_ids(edge_data)
        self.vertices = tuple(
            data.get_vertex() for data in vertex_data
        )
        self.edges = tuple(
            data.get_edge() for data in edge_data
        )
        self._construction_stages = construction_stages
        self._set_stages()
        self._num_processes = num_processes

    def _set_data_ids(self, data):
        for i, data in enumerate(data):
            data.id = i

    def _set_stages(self):
        self._stages = tuple(
            [] for i in range(len(self._construction_stages)+1)
        )
        for vertex in self.vertices:
            placed = False
            for i, stage in enumerate(self._construction_stages):
                if stage(vertex):
                    self._stages[i].append(vertex)
                    placed = True
                    break
            if not placed:
                self._stages[-1].append(vertex)


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

    def _place_building_blocks(self, mol, vertices, edges):
        self._implementation._place_building_blocks(molecule,)

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
