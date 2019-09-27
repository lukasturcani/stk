"""
Metal Complex
=============

#. :class:`.SquarePlanar`

"""


import logging
import numpy as np

from .topology_graph import TopologyGraph, VertexData, Vertex, EdgeData
from ...utilities import vector_angle

logger = logging.getLogger(__name__)


class _MetalComplexVertexData(VertexData):
    """
    Holds the data of a metal complex vertex.

    Attributes
    ----------
    id : :class:`int`
        The id of the vertex. Must match the index in
        :attr:`TopologyGraph.vertices`.

    position : :class:`numpy.ndarray`
        The position of the vertex.

    edges : :class:`list` of :class:`.EdgeData`
        The edges connected to the vertex.

    cell : :class:`numpy.ndarray`
        The unit cell in which the vertex is found.

    aligner_edge : :class:`int`
        The edge which is used to align the :class:`.BuildingBlock`
        placed on the vertex. The first :class:`.FunctionalGroup`
        in :attr:`.BuildingBlock.func_groups` is rotated such that
        it lies exactly on this :class:`.Edge`. Must be between
        ``0`` and the number of edges the vertex is connected to.

    """

    def __init__(self, x, y, z):
        """
        Initialize a :class:`_CageVertexData` instance.

        Parameters
        ----------
        x : :class:`float`
            The x coordinate.

        y : :class:`float`
            The y coordinate.

        z : :class:`float`
            The z coordinate.

        """

        self.aligner_edge = None
        super().__init__(x, y, z)

    def clone(self, clear_edges=False):
        """
        Return a clone.

        Parameters
        ----------
        clear_edges : :class:`bool`, optional
            ``True`` if the clone should not be connected to any edges.

        Returns
        -------
        :class:`_CageVertexData`
            The clone.

        """

        clone = super().clone(clear_edges)
        clone.aligner_edge = self.aligner_edge
        return clone

    def get_vertex(self):
        return _MetalComplexVertex(self)


class _MetalComplexVertex(Vertex):
    """
    Represents a vertex of a :class:`.MetalComplex`.

    Attributes
    ----------
    id : :class:`int`
        The id of the vertex. This should be its index in
        :attr:`TopologyGraph.vertices`.

    """

    def __init__(self, data):
        self._aligner_edge = data.aligner_edge
        super().__init__(data)

    def clone(self, clear_edges=False):
        """
        Return a clone.

        Parameters
        ----------
        clear_edges : :class:`bool`, optional
            ``True`` if the clone should not be connected to any edges.

        Returns
        -------
        :class:`Vertex`
            The clone.

        """

        clone = super().clone(clear_edges)
        clone._aligner_edge = self._aligner_edge
        return clone

    def get_aligner_edge(self):
        return self._aligner_edge

    def place_building_block(self, building_block, vertices, edges):
        """
        Place `building_block` on the :class:`.Vertex`.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block molecule which is to be placed on the
            vertex.

        vertices : :class:`tuple` of :class:`.Vertex`
            All vertices in the topology graph. The index of each
            vertex must match its :class:`~.Vertex.id`.

        edges : :class:`tuple` of :class:`.Edge`
            All edges in the topology graph. The index of each
            edge must match its :class:`~.Edge.id`.

        Returns
        -------
        :class:`numpy.nadarray`
            The position matrix of `building_block` after being
            placed.

        """

        bb_fg_names = list(set((
            i.fg_type.name for i in building_block.func_groups
        )))

        if 'metal' in bb_fg_names:
            return self._place_metal_atom(
                building_block=building_block,
                vertices=vertices,
                edges=edges
            )
        elif len(building_block.func_groups) == 1:
            return self._place_cap_building_block(
                building_block=building_block,
                vertices=vertices,
                edges=edges
            )
        elif len(building_block.func_groups) == 2:
            return self._place_linear_building_block(
                building_block=building_block,
                vertices=vertices,
                edges=edges
            )
        return self._place_nonlinear_building_block(
            building_block=building_block,
            vertices=vertices,
            edges=edges
        )

    def _place_metal_atom(
        self,
        building_block,
        vertices,
        edges
    ):
        """
        Place `building_block` on the :class:`.Vertex`.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block molecule which is to be placed on the
            vertex.

        vertices : :class:`tuple` of :class:`.Vertex`
            All vertices in the topology graph. The index of each
            vertex must match its :class:`~.Vertex.id`.

        edges : :class:`tuple` of :class:`.Edge`
            All edges in the topology graph. The index of each
            edge must match its :class:`~.Edge.id`.

        Returns
        -------
        :class:`numpy.nadarray`
            The position matrix of `building_block` after being
            placed.

        """
        building_block.set_centroid(position=self._position)
        return building_block.get_position_matrix()

    def _place_cap_building_block(
        self,
        building_block,
        vertices,
        edges
    ):
        """
        Place `building_block` on the :class:`.Vertex`.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block molecule which is to be placed on the
            vertex.

        vertices : :class:`tuple` of :class:`.Vertex`
            All vertices in the topology graph. The index of each
            vertex must match its :class:`~.Vertex.id`.

        edges : :class:`tuple` of :class:`.Edge`
            All edges in the topology graph. The index of each
            edge must match its :class:`~.Edge.id`.

        Returns
        -------
        :class:`numpy.nadarray`
            The position matrix of `building_block` after being
            placed.

        """
        building_block.set_centroid(position=self._position)
        fg_centroid = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_bonder_ids()
        )
        start = fg_centroid - self._position
        aligner_edge = edges[self._edge_ids[self._aligner_edge]]
        edge_coord = aligner_edge.get_position()
        target = edge_coord - self._position
        building_block.apply_rotation_between_vectors(
            start=start,
            target=target,
            origin=building_block.get_centroid()
        )

        return building_block.get_position_matrix()

    def _place_linear_building_block(
        self,
        building_block,
        vertices,
        edges
    ):
        """
        Place `building_block` on the :class:`.Vertex`.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block molecule which is to be placed on the
            vertex.

        vertices : :class:`tuple` of :class:`.Vertex`
            All vertices in the topology graph. The index of each
            vertex must match its :class:`~.Vertex.id`.

        edges : :class:`tuple` of :class:`.Edge`
            All edges in the topology graph. The index of each
            edge must match its :class:`~.Edge.id`.

        Returns
        -------
        :class:`numpy.nadarray`
            The position matrix of `building_block` after being
            placed.

        """

        building_block.set_centroid(
            position=self._position,
            atom_ids=building_block.get_bonder_ids()
        )
        fg_centroid = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_bonder_ids()
        )
        start = fg_centroid - self._position
        aligner_edge = edges[self._edge_ids[self._aligner_edge]]
        edge_coord = aligner_edge.get_position()
        connected_edges = tuple(edges[id_] for id_ in self._edge_ids)
        target = edge_coord - self._get_edge_centroid(
            centroid_edges=connected_edges,
            vertices=vertices
        )
        building_block.apply_rotation_between_vectors(
            start=start,
            target=target,
            origin=self._position
        )
        start = building_block.get_centroid_centroid_direction_vector()
        e0_coord = edges[self._edge_ids[0]].get_position()
        e1_coord = edges[self._edge_ids[1]].get_position()
        building_block.apply_rotation_to_minimize_angle(
            start=start,
            target=self._position,
            axis=e0_coord-e1_coord,
            origin=self._position,
        )
        return building_block.get_position_matrix()

    def _place_nonlinear_building_block(
        self,
        building_block,
        vertices,
        edges
    ):
        """
        Place `building_block` on the :class:`.Vertex`.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block molecule which is to be placed on the
            vertex.

        vertices : :class:`tuple` of :class:`.Vertex`
            All vertices in the topology graph. The index of each
            vertex must match its :class:`~.Vertex.id`.

        edges : :class:`tuple` of :class:`.Edge`
            All edges in the topology graph. The index of each
            edge must match its :class:`~.Edge.id`.

        Returns
        -------
        :class:`numpy.nadarray`
            The position matrix of `building_block` after being
            placed.

        """

        building_block.set_centroid(
            position=self._position,
            atom_ids=building_block.get_bonder_ids()
        )
        connected_edges = tuple(edges[id_] for id_ in self._edge_ids)
        edge_normal = self._get_edge_plane_normal(
            reference=self._get_edge_centroid(
                centroid_edges=connected_edges,
                vertices=vertices
            ),
            plane_edges=connected_edges,
            vertices=vertices
        )
        building_block.apply_rotation_between_vectors(
            start=building_block.get_bonder_plane_normal(),
            target=edge_normal,
            origin=self._position
        )
        fg_bonder_centroid = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_bonder_ids()
        )
        start = fg_bonder_centroid - self._position
        aligner_edge = edges[self._edge_ids[self._aligner_edge]]
        edge_coord = aligner_edge.get_position()
        target = edge_coord - self._get_edge_centroid(
            centroid_edges=connected_edges,
            vertices=vertices
        )
        building_block.apply_rotation_to_minimize_angle(
            start=start,
            target=target,
            axis=edge_normal,
            origin=self._position
        )
        return building_block.get_position_matrix()

    def assign_func_groups_to_edges(
        self,
        building_block,
        vertices,
        edges
    ):
        """
        Assign functional groups to edges.

        Each :class:`.FunctionalGroup` of the `building_block` needs
        to be associated with one of the :class:`.Edge` instances in
        :attr:`edges`.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is needs to have
            functional groups assigned to edges.

        vertices : :class:`tuple` of :class:`.Vertex`
            All vertices in the topology graph. The index of each
            vertex must match its :class:`~.Vertex.id`.

        edges : :class:`tuple` of :class:`.Edge`
            All edges in the topology graph. The index of each
            edge must match its :class:`~.Edge.id`.

        Returns
        -------
        :class:`dict`
            A mapping from the id of a functional group in
            `building_block` to the id of the edge in :attr:`edges` it
            is assigned to.

        """

        bb_fg_names = list(set((
            i.fg_type.name for i in building_block.func_groups
        )))
        if 'metal' in bb_fg_names:
            return self._assign_func_groups_to_metal_atom(
                building_block=building_block,
                vertices=vertices,
                edges=edges
            )
        elif len(building_block.func_groups) == 1:
            return self._assign_func_groups_to_cap_edges(
                building_block=building_block,
                vertices=vertices,
                edges=edges
            )
        elif len(building_block.func_groups) == 2:
            return self._assign_func_groups_to_linear_edges(
                building_block=building_block,
                vertices=vertices,
                edges=edges
            )
        return self._assign_func_groups_to_nonlinear_edges(
            building_block=building_block,
            vertices=vertices,
            edges=edges
        )

    def after_assign_func_groups_to_edges(
        self,
        building_block,
        func_groups,
        vertices,
        edges
    ):
        """
        Perform operations after functional groups have been assigned.

        This method is always executed serially. It is often useful
        when data needs to be transferred between vertices, which
        have been processed independently, in parallel.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is needs to have
            functional groups assigned to edges.

        func_groups : :class:`tuple` of :class:`.FunctionalGroup`
            The functional group clones added to the constructed
            molecule.

        vertices : :class:`tuple` of :class:`.Vertex`
            All vertices in the topology graph. The index of each
            vertex must match its :class:`~.Vertex.id`.

        edges : :class:`tuple` of :class:`.Edge`
            All edges in the topology graph. The index of each
            edge must match its :class:`~.Edge.id`.

        Returns
        -------
        None : :class:`NoneType`

        """

        bb_fgs = set(func_groups)
        for edge_id in self._edge_ids:
            for func_group in edges[edge_id].get_func_groups():
                if func_group not in bb_fgs:
                    continue

                for vertex_id in edges[edge_id].get_vertex_ids():
                    if vertex_id == self.id:
                        continue

        return super().after_assign_func_groups_to_edges(
            building_block=building_block,
            func_groups=func_groups,
            vertices=vertices,
            edges=edges
        )

    def _assign_func_groups_to_metal_atom(
        self,
        building_block,
        vertices,
        edges
    ):
        return {
            fg_id: edge_id for fg_id, edge_id in enumerate(sorted(
                self._edge_ids,
                key=self._get_fg0_distance(building_block, edges)
            ))
        }

    def _assign_func_groups_to_cap_edges(
        self,
        building_block,
        vertices,
        edges
    ):
        return {
            fg_id: edge_id for fg_id, edge_id in enumerate(sorted(
                self._edge_ids,
                key=self._get_fg0_distance(building_block, edges)
            ))
        }

    def _assign_func_groups_to_linear_edges(
        self,
        building_block,
        vertices,
        edges
    ):
        return {
            fg_id: edge_id for fg_id, edge_id in enumerate(sorted(
                self._edge_ids,
                key=self._get_fg0_distance(building_block, edges)
            ))
        }

    def _get_fg0_distance(self, building_block, edges):
        fg_coord = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_bonder_ids()
        )

        def distance(edge_id):
            displacement = edges[edge_id].get_position() - fg_coord
            return np.linalg.norm(displacement)

        return distance

    def _assign_func_groups_to_nonlinear_edges(
        self,
        building_block,
        vertices,
        edges
    ):
        # The idea is to order the functional groups in building_block
        # by their angle from func_groups[0] and the bonder centroid,
        #  going in the clockwise direction.
        #
        # The edges are also ordered by their angle from aligner_edge
        # and the edge centroid going in the clockwise direction.
        #
        # Once the fgs and edges are ordered, zip and assign them.

        fg0_coord = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_bonder_ids()
        )
        bonder_centroid = building_block.get_centroid(
            atom_ids=building_block.get_bonder_ids()
        )
        fg0_direction = fg0_coord-bonder_centroid
        axis = np.cross(
            fg0_direction,
            building_block.get_bonder_plane_normal()
        )
        func_groups = sorted(
            range(len(building_block.func_groups)),
            key=self._get_func_group_angle(
                building_block=building_block,
                fg0_direction=fg0_direction,
                bonder_centroid=bonder_centroid,
                axis=axis
            )
        )
        assignments = {}
        edge_ids = sorted(
            self._edge_ids,
            key=self._get_edge_angle(axis, vertices, edges)
        )
        for edge_id, fg_id in zip(edge_ids, func_groups):
            assignments[fg_id] = edge_id
        return assignments

    @staticmethod
    def _get_func_group_angle(
        building_block,
        fg0_direction,
        bonder_centroid,
        axis
    ):

        def angle(fg_id):
            func_group = building_block.func_groups[fg_id]
            coord = building_block.get_centroid(
                atom_ids=func_group.get_bonder_ids()
            )
            fg_direction = coord-bonder_centroid
            theta = vector_angle(fg0_direction, fg_direction)

            projection = fg_direction @ axis
            if theta > 0 and projection < 0:
                return 2*np.pi - theta
            return theta

        return angle

    def _get_edge_angle(self, axis, vertices, edges):
        aligner_edge = edges[self._edge_ids[self._aligner_edge]]
        aligner_edge_coord = aligner_edge.get_position()
        connected_edges = tuple(edges[id_] for id_ in self._edge_ids)
        edge_centroid = self._get_edge_centroid(
            centroid_edges=connected_edges,
            vertices=vertices
        )
        # This axis is used to figure out the clockwise direction.
        aligner_edge_direction = aligner_edge_coord - edge_centroid

        def angle(edge_id):
            coord = edges[edge_id].get_position()
            edge_direction = coord - edge_centroid
            theta = vector_angle(
                vector1=edge_direction,
                vector2=aligner_edge_direction
            )

            projection = edge_direction @ axis
            if theta > 0 and projection < 0:
                return 2*np.pi - theta
            return theta

        return angle

    def __str__(self):
        return (
            f'Vertex(id={self.id}, '
            f'position={self._position.tolist()}, '
            f'aligner_edge={self._aligner_edge})'
        )


class MetalComplex(TopologyGraph):
    """
    Represents single-molecule metal complex topology graphs.

    MetalComplex topologies are added by creating a subclass which
    defines the :attr:`vertex_data` and :attr:`edge_data` class
    attributes.

    This class is modelled after :class:`Cage` and its subclasses.

    Attributes
    ----------
    vertex_data : :class:`tuple` of :class:`.VertexData`
        A class attribute. Holds the data of the vertices which make up
        the topology graph.

    edge_data : :class:`tuple` of :class:`.EdgeData`
        A class attribute. Holds the data of the edges which make up
        the topology graph.

    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    Examples
    --------


    """

    def __init_subclass__(cls, **kwargs):
        for i, vertex in enumerate(cls.vertex_data):
            vertex.id = i
        for i, edge in enumerate(cls.edge_data):
            edge.id = i
        return super().__init_subclass__(**kwargs)

    def __init__(self, vertex_alignments=None,
                 unsatured_vertices=None, num_processes=1):
        """
        Initialize a :class:`.MetalComplex`.

        Parameters
        ----------
        vertex_alignments : :class:`dict`, optional
            A mapping from the :attr:`.Vertex.id` of a :class:`.Vertex`
            :attr:`vertices` to an :class:`.Edge` connected to it.
            The :class:`.Edge` is used to align the first
            :class:`.FunctionalGroup` of a :class:`.BuildingBlock`
            placed on that vertex. Only vertices which need to have
            their default edge changed need to be present in the
            :class:`dict`. If ``None`` then the default edge is used
            for each vertex. Changing which :class:`.Edge` is used will
            mean that the topology graph represents different
            structural isomers. The edge is refered to by a number
            between ``0`` (inclusive) and the number of edges the
            vertex is connected to (exclusive).

        unsaturated_vertices : :class:`list` of :class:`int`, optional
            A list of the unsaturated sites on the metal complexes to
            be built. The integers correspond to vertex ids.

        num_processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        """

        # Metal complexes can have unsaturated sites.
        # Need to remove information about the sites that will not
        # react from stage, self.vertices and self.edges.
        if unsatured_vertices is not None:
            self.old_vertex_data = self.vertex_data
            self.old_edge_data = self.edge_data
            self.vertex_data = tuple(
                i for i in self.old_vertex_data
                if i.id not in unsatured_vertices
            )
            used_edges = [
                i for i in self.old_edge_data
                if set(i.vertices).issubset(set(self.vertex_data))
            ]
            self.edge_data = tuple(i for i in used_edges)

        if vertex_alignments is None:
            vertex_alignments = {}

        vertex_data = {
            data: data.clone(True) for data in self.vertex_data
        }
        for vertex in vertex_data.values():
            vertex.aligner_edge = vertex_alignments.get(vertex.id, 0)
        edge_data = tuple(
            edge.clone(vertex_data)
            for edge in self.edge_data
        )
        vertex_types = sorted(
            {len(v.edges) for v in vertex_data},
            reverse=True
        )
        super().__init__(
            vertex_data=vertex_data.values(),
            edge_data=edge_data,
            construction_stages=tuple(
                lambda vertex, vertex_type=vt:
                    vertex.get_num_edges() == vertex_type
                for vt in vertex_types
            ),
            num_processes=num_processes
        )

    def _get_scale(self, mol):
        """
        Get the scale used for the positions of :attr:`vertices`.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The molecule being constructed.

        Returns
        -------
        :class:`float` or :class:`list` of :class:`float`
            The value by which the position of each :class:`Vertex` is
            scaled. Can be a single number if all axes are scaled by
            the same amount or a :class:`list` of three numbers if
            each axis is scaled by a different value.

        """
        organic_bbs = [
            i for i in mol.building_block_vertices
            if 'metal' not in list(set((
                j.fg_type.name for j in i.func_groups
            )))
        ]
        if organic_bbs:
            return max(bb.get_maximum_diameter() for bb in organic_bbs)
        else:
            # No organic building blocks.
            return 1

    def __repr__(self):
        vertex_alignments = ', '.join(
            f'{v.id}: {v.get_aligner_edge()}'
            for v in self.vertices
        )
        return (
            f'metal_complex.{self.__class__.__name__}('
            f'vertex_alignments={{{vertex_alignments}}})'
        )


class SquarePlanarMonodentate(MetalComplex):
    """
    Represents a square planar metal complex topology graph.

    See :class:`.MetalComplex` for more details and examples.

    Attributes
    ----------
    vertex_data : :class:`tuple` of :class:`.VertexData`
        A class attribute. Holds the data of the vertices which make up
        the topology graph.

    edge_data : :class:`tuple` of :class:`.EdgeData`
        A class attribute. Holds the data of the edges which make up
        the topology graph.

    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    vertex_data = (
        _MetalComplexVertexData(0, 0, 0),
        _MetalComplexVertexData(0, 1, 0),
        _MetalComplexVertexData(0, 0, 1),
        _MetalComplexVertexData(0, -1, 0),
        _MetalComplexVertexData(0, 0, -1),
    )

    edge_data = (
        EdgeData(
            vertex_data[0],
            vertex_data[1],
            position=[0, 0.2, 0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[2],
            position=[0, 0, 0.2]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[3],
            position=[0, -0.2, 0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[4],
            position=[0, 0, 0.2]
        ),
    )


class SquarePlanarBidentate(MetalComplex):
    """
    Represents a square planar metal complex topology graph.

    See :class:`.MetalComplex` for more details and examples.

    Attributes
    ----------
    vertex_data : :class:`tuple` of :class:`.VertexData`
        A class attribute. Holds the data of the vertices which make up
        the topology graph.

    edge_data : :class:`tuple` of :class:`.EdgeData`
        A class attribute. Holds the data of the edges which make up
        the topology graph.

    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    vertex_data = (
        _MetalComplexVertexData(0, 0, 0),
        _MetalComplexVertexData(0, 1, 0),
        _MetalComplexVertexData(0, -1, 0)
    )

    edge_data = (
        EdgeData(
            vertex_data[0],
            vertex_data[1],
            position=[0.2, 0.2, 0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[1],
            position=[-0.2, 0.2, 0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[2],
            position=[0.2, -0.2, 0]
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[2],
            position=[-0.2, -0.2, 0]
        )
    )
