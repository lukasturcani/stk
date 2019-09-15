"""
Metal Complex
=============

#. :class:`.SquarePlanar`

"""


import logging
import numpy as np

from .topology_graph import TopologyGraph, Vertex, Edge
from ...utilities import vector_angle

logger = logging.getLogger(__name__)


class _MetalComplexVertex(Vertex):
    """
    Represents a vertex of a :class:`.MetalComplex`.

    Attributes
    ----------
    id : :class:`int`
        The id of the vertex. This should be its index in
        :attr:`TopologyGraph.vertices`.

    edges : :class:`list` of :class:`.Edge`
        The edges the :class:`Vertex` is connected to.

    aligner_edge : :class:`.Edge`
        The :class:`.Edge` in :attr:`edges`, which is used to align the
        :class:`.BuildingBlock` placed on the vertex. The first
        :class:`.FunctionalGroup` in :attr:`.BuildingBlock.func_groups`
        is rotated such that it lies exactly on this :class:`.Edge`.

    """

    def __init__(self, x, y, z, use_bonder_placement=True):
        """
        Initialize a :class:`_MetalComplexVertex`.

        Parameters
        ----------
        x : :class:`float`
            The x coordinate.

        y : :class:`float`
            The y coordinate.

        z : :class:`float`
            The z coordinate.

        use_bonder_placement : :class:`bool`, optional
            If ``True``the position of the vertex will be updated such
            that it is in the middle of the neighboring bonder
            centroids, rather than in the middle of the neighboring
            vertices.

        """

        # _neighbor_positions holds the bonder centroids of functional
        # groups on neighbor vertices connected to this vertex.
        self._neighbor_positions = []
        self.aligner_edge = None
        self._use_bonder_placement = use_bonder_placement
        super().__init__(x, y, z)

    @classmethod
    def init_at_center(cls, *vertices):
        """
        Initialize at the center of `vertices`.

        Parameters
        ----------
        vertices : :class:`.Vertex`
            Vertices at whose center this vertex should be initialized.

        Returns
        -------
        :class:`.Vertex`
            The vertex.

        """

        center = sum(vertex.get_position() for vertex in vertices)
        center /= len(vertices)
        return cls(*center)

    def clone(self, clear_edges=False):
        """
        Create a clone of the instance.

        Parameters
        ----------
        clear_edges : :class:`bool`, optional
            If ``True`` the :attr:`edges` attribute of the clone will
            be empty.

        Returns
        -------
        :class:`Vertex`
            A clone with the same position but not connected to any
            :class:`.Edge` objects.

        """

        clone = super().clone(clear_edges)
        if self.aligner_edge is None:
            clone.aligner_edge = None
        else:
            clone.aligner_edge = self.aligner_edge.clone(
                add_to_vertices=False
            )
        clone._use_bonder_placement = self._use_bonder_placement
        clone._neighbor_positions = list(self._neighbor_positions)
        return clone

    def apply_scale(self, scale):
        """
        Scale the position by `scale`.

        Parameters
        ----------
        scale : :class:`float` or :class:`list`of :class:`float`
            The value by which the position of the :class:`Vertex` is
            scaled. Can be a single number if all axes are scaled by
            the same amount or a :class:`list` of three numbers if
            each axis is scaled by a different value.

        Returns
        -------
        :class:`Vertex`
            The vertex is returned.

        """

        self._position *= scale
        self.aligner_edge.apply_scale(scale)
        return self

    def place_building_block(self, building_block):
        """
        Place `building_block` on the :class:`.Vertex`.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block molecule which is to be placed on the
            vertex.

        Returns
        -------
        :class:`numpy.nadarray`
            The position matrix of `building_block` after being
            placed.

        """

        if (
            self._use_bonder_placement
            and len(self._neighbor_positions) == len(self.edges)
        ):
            self._update_position()

        bb_fg_names = list(set((
            i.fg_type.name for i in building_block.func_groups
        )))

        if bb_fg_names[0] == 'metal':
            return self._place_metal_atom(building_block)
        elif len(building_block.func_groups) == 1:
            return self._place_cap_building_block(building_block)
        elif len(building_block.func_groups) == 2:
            return self._place_linear_building_block(building_block)
        return self._place_nonlinear_building_block(building_block)

    def _update_position(self):
        self._position = np.divide(
            np.sum(self._neighbor_positions, axis=0),
            len(self._neighbor_positions)
        )

    def _place_metal_atom(self, building_block):
        """
        Place `building_block` on the :class:`.Vertex`.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block molecule which is to be placed on the
            vertex.

        Returns
        -------
        :class:`numpy.nadarray`
            The position matrix of `building_block` after being
            placed.

        """

        building_block.set_centroid(
            position=self._position
        )
        print('placed', building_block.get_position_matrix())
        return building_block.get_position_matrix()

    def _place_cap_building_block(self, building_block):
        """
        Place `building_block` on the :class:`.Vertex`.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block molecule which is to be placed on the
            vertex.

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
        print('ff', fg_centroid)
        start = fg_centroid - self._position
        print(self._position)
        print('s', start)
        print(self.aligner_edge)
        edge_coord = self.aligner_edge.get_position()
        print('e', edge_coord)
        print('ge', self._get_edge_centroid())
        target = edge_coord
        print('t', target)
        building_block.apply_rotation_between_vectors(
            start=start,
            target=target,
            origin=self._position
        )
        input()
        return building_block.get_position_matrix()

    def _place_linear_building_block(self, building_block):
        """
        Place `building_block` on the :class:`.Vertex`.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block molecule which is to be placed on the
            vertex.

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
        edge_coord = self.aligner_edge.get_position()
        target = edge_coord - self._get_edge_centroid()
        building_block.apply_rotation_between_vectors(
            start=start,
            target=target,
            origin=self._position
        )
        start = building_block.get_centroid_centroid_direction_vector()
        e0_coord = self.edges[0].get_position()
        e1_coord = self.edges[1].get_position()
        building_block.apply_rotation_to_minimize_angle(
            start=start,
            target=self._position,
            axis=e0_coord-e1_coord,
            origin=self._position,
        )
        return building_block.get_position_matrix()

    def _place_nonlinear_building_block(self, building_block):
        """
        Place `building_block` on the :class:`.Vertex`.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block molecule which is to be placed on the
            vertex.

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
        edge_normal = self._get_edge_plane_normal(
            reference=self._get_edge_centroid()
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
        edge_coord = self.aligner_edge.get_position()
        target = edge_coord - self._get_edge_centroid()
        building_block.apply_rotation_to_minimize_angle(
            start=start,
            target=target,
            axis=edge_normal,
            origin=self._position
        )
        return building_block.get_position_matrix()

    def assign_func_groups_to_edges(self, building_block):
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
        if bb_fg_names[0] == 'metal':
            return self._assign_func_groups_to_metal_atom(
                building_block=building_block
            )
        elif len(building_block.func_groups) == 1:
            return self._assign_func_groups_to_cap_edges(
                building_block=building_block
            )
        elif len(building_block.func_groups) == 2:
            return self._assign_func_groups_to_linear_edges(
                building_block=building_block
            )
        return self._assign_func_groups_to_nonlinear_edges(
            building_block=building_block
        )

    def after_assign_func_groups_to_edges(
        self,
        building_block,
        func_groups
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

        Returns
        -------
        None : :class:`NoneType`

        """

        bb_fgs = set(func_groups)
        for edge in self.edges:
            for func_group in edge.get_func_groups():
                if func_group not in bb_fgs:
                    continue

                bonder_position = self._get_molecule_centroid(
                    atom_ids=func_group.get_bonder_ids()
                )
                for vertex in edge.vertices:
                    if vertex is self:
                        continue

                    vertex._neighbor_positions.append(bonder_position)

        return super().after_assign_func_groups_to_edges(
            building_block=building_block,
            func_groups=func_groups
        )

    def _assign_func_groups_to_metal_atom(self, building_block):
        return {
            fg_id: e.id for fg_id, e in enumerate(sorted(
                self.edges,
                key=self._get_fg0_distance(building_block)
            ))
        }

    def _assign_func_groups_to_cap_edges(self, building_block):
        return {
            fg_id: e.id for fg_id, e in enumerate(sorted(
                self.edges,
                key=self._get_fg0_distance(building_block)
            ))
        }

    def _assign_func_groups_to_linear_edges(self, building_block):
        return {
            fg_id: e.id for fg_id, e in enumerate(sorted(
                self.edges,
                key=self._get_fg0_distance(building_block)
            ))
        }

    def _get_fg0_distance(self, building_block):
        fg_coord = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_bonder_ids()
        )

        def distance(edge):
            displacement = edge.get_position() - fg_coord
            return np.linalg.norm(displacement)

        return distance

    def _assign_func_groups_to_nonlinear_edges(self, building_block):
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
        edges = sorted(self.edges, key=self._get_edge_angle(axis))
        for edge, fg_id in zip(edges, func_groups):
            assignments[fg_id] = edge.id
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

    def _get_edge_angle(self, axis):

        aligner_edge_coord = self.aligner_edge.get_position()
        edge_centroid = self._get_edge_centroid()
        # This axis is used to figure out the clockwise direction.
        aligner_edge_direction = aligner_edge_coord - edge_centroid

        def angle(edge):
            coord = edge.get_position()
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


class MetalComplex(TopologyGraph):
    """
    Represents single-molecule metal complex topology graphs.

    MetalComplex topologies are added by creating a subclass which
    defines the :attr:`vertices` and :attr:`edges` of the topology
    as class attributes.

    This class is modelled after :class:`Cage` and its subclasses.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    Examples
    --------


    """

    def __init_subclass__(cls, **kwargs):
        for i, vertex in enumerate(cls.vertices):
            vertex.id = i
        return super().__init_subclass__(**kwargs)

    def __init__(self, vertex_alignments=None,
                 unsatured_vertices=None, num_processes=1):
        """
        Initialize a :class:`.MetalComplex`.

        Parameters
        ----------
        vertex_alignments : :class:`dict`, optional
            A mapping from a :class:`.Vertex` in :attr:`vertices`
            to an :class:`.Edge` connected to it. The :class:`.Edge` is
            used to align the first :class:`.FunctionalGroup` of a
            :class:`.BuildingBlock` placed on that vertex. Only
            vertices which need to have their default edge changed need
            to be present in the :class:`dict`. If ``None`` then the
            first :class:`.Edge` in :class:`.Vertex.edges` is for each
            vertex is used. Changing which :class:`.Edge` is used will
            mean that the topology graph represents different
            structural isomers.

            The vertices and edges can also be referred to by their
            indices.
        num_processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        """
        if vertex_alignments is None:
            vertex_alignments = {}

        # Convert ints to Vertex and Edge instances.
        _vertex_alignments = {}
        for v, e in vertex_alignments.items():
            v = self.vertices[v] if isinstance(v, int) else v
            e = v.edges[e] if isinstance(e, int) else e
            _vertex_alignments[v] = e
        vertex_alignments = _vertex_alignments

        vertex_clones = {}
        for vertex in self.vertices:
            clone = vertex.clone(clear_edges=True)
            clone.aligner_edge = vertex_alignments.get(
                vertex,
                vertex.edges[0]
            )
            vertex_clones[vertex] = clone

        edge_clones = {}
        for edge in self.edges:
            edge_clones[edge] = edge.clone(vertex_clones)

        vertices = tuple(vertex_clones.values())
        for vertex in vertices:
            vertex.aligner_edge = edge_clones[vertex.aligner_edge]

        vertex_types = sorted(
            set(len(v.edges) for v in self.vertices),
            reverse=True
        )
        super().__init__(
            vertices=vertices,
            edges=tuple(edge_clones.values()),
            construction_stages=tuple(
                lambda vertex, vertex_type=vt:
                    len(vertex.edges) == vertex_type
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

        return max(
            bb.get_maximum_diameter()
            for bb in mol.building_block_vertices
        )

    def _clean_up(self, mol):
        print(
            'here you should determine how many FGs remain, and assign'
        )
        return super()._clean_up(mol)

    def __repr__(self):
        vertex_alignments = ', '.join(
            f'{v.id}: {v.edges.index(v.aligner_edge)}'
            for v in self.vertices
        )
        return (
            f'MetalComplex.{self.__class__.__name__}('
            f'vertex_alignments={{{vertex_alignments}}})'
        )


class SquarePlanar(MetalComplex):
    """
    Represents a square planar metal complex topology graph.

    See :class:`.MetalComplex` for more details and examples.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    vertices = (
        _MetalComplexVertex(0, 0, 0),
        _MetalComplexVertex(0, 1, 0),
        _MetalComplexVertex(0, 0, 1),
        _MetalComplexVertex(0, -1, 0),
        _MetalComplexVertex(0, 0, -1),
    )

    edges = (
        Edge(vertices[0], vertices[1], position=[0, 0.5, 0]),
        Edge(vertices[0], vertices[2], position=[0, 0, 0.5]),
        Edge(vertices[0], vertices[3], position=[0, -0.5, 0]),
        Edge(vertices[0], vertices[4], position=[0, 0, 0.5]),
    )
