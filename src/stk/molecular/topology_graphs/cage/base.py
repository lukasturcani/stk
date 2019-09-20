"""
Organic Cage
============

For usage examples see :class:`.Cage`.

#. :class:`.TwoPlusThree`
#. :class:`.FourPlusSix`
#. :class:`.FourPlusSix2`
#. :class:`.SixPlusNine`
#. :class:`.EightPlusTwelve`
#. :class:`.TwentyPlusThirty`
#. :class:`.TwoPlusFour`
#. :class:`.ThreePlusSix`
#. :class:`.FourPlusEight`
#. :class:`.FivePlusTen`
#. :class:`.SixPlusTwelve`
#. :class:`.EightPlusSixteen`
#. :class:`.TenPlusTwenty`
#. :class:`.SixPlusEight`
#. :class:`.OnePlusOne`
#. :class:`.TwoPlusTwo`
#. :class:`.FourPlusFour`
#. :class:`.TwelvePlusThirty`

"""


import numpy as np

from ..topology_graph import TopologyGraph, VertexData, Vertex
from ....utilities import vector_angle


class _CageVertexData(VertexData):
    """
    Holds the data of a cage vertex.

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

    use_bonder_placement : :class:`bool`, optional
        If ``True``the position of the vertex will be updated such
        that it is in the middle of the neighboring bonder
        centroids, rather than in the middle of the neighboring
        vertices.

    """

    def __init__(self, x, y, z, use_bonder_placement=True):
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

        use_bonder_placement : :class:`bool`, optional
            If ``True``the position of the vertex will be updated such
            that it is in the middle of the neighboring bonder
            centroids, rather than in the middle of the neighboring
            vertices.

        """

        self.aligner_edge = None
        self.use_bonder_placement = use_bonder_placement
        super().__init__(x, y, z)

    @classmethod
    def init_at_center(cls, *vertex_data):
        obj = super().init_at_center(*vertex_data)
        obj.aligner_edge = None
        obj.use_bonder_placement = True
        return obj

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
        clone.use_bonder_placement = self.use_bonder_placement
        return clone

    def get_vertex(self):
        return _CageVertex(self)


class _CageVertex(Vertex):
    """
    Represents a vertex of a :class:`.Cage`.

    Attributes
    ----------
    id : :class:`int`
        The id of the vertex. This should be its index in
        :attr:`TopologyGraph.vertices`.

    """

    def __init__(self, data):
        # _neighbor_positions holds the bonder centroids of functional
        # groups on neighbor vertices connected to this vertex.
        self._neighbor_positions = []
        self._use_bonder_placement = data.use_bonder_placement
        # The edge which is used to align the :class:`.BuildingBlock`
        # placed on the vertex. The first :class:`.FunctionalGroup`
        # in :attr:`.BuildingBlock.func_groups` is rotated such that
        # it lies exactly on this :class:`.Edge`. Must be between
        # ``0`` and the number of edges the vertex is connected to.
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
        clone._use_bonder_placement = self._use_bonder_placement
        clone._neighbor_positions = list(self._neighbor_positions)
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

        if (
            self._use_bonder_placement
            and len(self._neighbor_positions) == len(self._edge_ids)
        ):
            self._update_position()

        if len(building_block.func_groups) == 2:
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

    def _update_position(self):
        self._position = np.divide(
            np.sum(self._neighbor_positions, axis=0),
            len(self._neighbor_positions)
        )

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

        if len(building_block.func_groups) == 2:
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

                bonder_position = self._get_molecule_centroid(
                    atom_ids=func_group.get_bonder_ids()
                )
                for vertex_id in edges[edge_id].get_vertex_ids():
                    if vertex_id == self.id:
                        continue

                    vertices[vertex_id]._neighbor_positions.append(
                        bonder_position
                    )

        return super().after_assign_func_groups_to_edges(
            building_block=building_block,
            func_groups=func_groups,
            vertices=vertices,
            edges=edges
        )

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


class Cage(TopologyGraph):
    """
    Represents a cage topology graph.

    Cage topologies are added by creating a subclass which defines the
    :attr:`vertex_data` and :attr:`edge_data` class attributes.

    A :class:`Cage` subclass will add the attributes
    :attr:`num_windows` and :attr:`num_window_types` to each
    :class:`.ConstructedMolecule`.

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
    :class:`Cage` instances can be made without supplying
    additional arguments (using :class:`.FourPlusSix` as an example)

    .. code-block:: python

        import stk

        bb1 = stk.BuildingBlock('NCCN', ['amine'])
        bb2 = stk.BuildingBlock('O=CC(C=O)C=O', ['aldehyde'])
        cage1 = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topology_graph=stk.cage.FourPlusSix()
        )

    Different structural isomers of cages can be made by using the
    `vertex_alignments` optional parameter

    .. code-block:: python

        tetrahedron = stk.cage.FourPlusSix(
            vertex_alignments={0: 1, 1: 1, 2: 2}
        )
        cage2 = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topology_graph=tetrahedron
        )

    The parameter maps the :attr:`~.Vertex.id` of a vertex to a number
    between 0 (inclusive) and the number of edges the vertex is
    connected to (exclusive). So a vertex connected to three edges
    can be mapped to ``0``, ``1`` or ``2``.

    By changing which edge each vertex is aligned with, a different
    structural isomer of the cage can be formed.

    You can also build cages with multiple building blocks, but you
    have to assign each building block to a vertex with
    `building_block_vertices`.

    bb1 = stk.BuildingBlock('O=CC(C=O)C=O', ['aldehyde'])
    bb2 = stk.BuildingBlock('O=CC(Cl)(C=O)C=O', ['aldehyde'])
    bb3 = stk.BuildingBlock('NCCN', ['amine'])
    bb4 = stk.BuildingBlock('NCC(Cl)N', ['amine'])
    bb5 = stk.BuildingBlock('NCCCCN', ['amine'])

    tetrahedron = stk.cage.FourPlusSix()
    cage = stk.ConstructedMolecule(
        building_blocks=[bb1, bb2, bb3, bb4, bb5],
        topology_graph=tetrahedron,
        building_block_vertices={
            bb1: tetrahedron.vertices[:2],
            bb2: tetrahedron.vertices[2:4],
            bb3: tetrahedron.vertices[4:5],
            bb4: tetrahedron.vertices[5:6],
            bb5: tetrahedron.vertices[6:]
        }
    )

    """

    def __init_subclass__(cls, **kwargs):
        for i, vertex in enumerate(cls.vertex_data):
            vertex.id = i
        for i, edge in enumerate(cls.edge_data):
            edge.id = i
        return super().__init_subclass__(**kwargs)

    def __init__(self, vertex_alignments=None, num_processes=1):
        """
        Initialize a :class:`.Cage`.

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

        num_processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        """

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

    def assign_building_blocks_to_vertices(self, building_blocks):
        """
        Assign `building_blocks` to :attr:`vertices`.

        Parameters
        ----------
        building_blocks : :class:`list` of :class:`.Molecule`
            The :class:`.BuildingBlock` and
            :class:`ConstructedMolecule` instances which
            represent the building block molecules used for
            construction. Only one instance is present per building
            block molecule, even if multiples of that building block
            join up to form the :class:`ConstructedMolecule`.

        Returns
        -------
        :class:`dict`
            Maps the `building_blocks`, to the
            :class:`~.topologies.base.Vertex` objects in
            :attr:`vertices` they are placed on during construction.
            The :class:`dict` has the form

            .. code-block:: python

                building_block_vertices = {
                    BuildingBlock(...): [Vertex(...), Vertex(...)],
                    BuildingBlock(...): [
                        Vertex(...),
                        Vertex(...),
                        Vertex(...),
                    ]
                    ConstructedMolecule(...): [Vertex(...)]
                }

        Raises
        ------
        :class:`ValueError`
            If there is more than one building with a given number
            of functional groups.

        """

        bb_by_degree = {}
        for bb in building_blocks:
            num_fgs = len(bb.func_groups)
            if num_fgs in bb_by_degree:
                raise ValueError(
                    'If there are multiple building blocks with the '
                    'same number of functional groups, '
                    'building_block_vertices must be set explicitly.'
                )
            bb_by_degree[num_fgs] = bb

        building_block_vertices = {}
        for vertex in self.vertices:
            bb = bb_by_degree[vertex.get_num_edges()]
            building_block_vertices[bb] = (
                building_block_vertices.get(bb, [])
            )
            building_block_vertices[bb].append(vertex)
        return building_block_vertices

    def _prepare(self, mol):
        """
        Do preprocessing on `mol` before construction.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The molecule being constructed.

        Returns
        -------
        None : :class:`NoneType`

        """

        # Order the building blocks by number of functional groups
        # so that building blocks with more functional groups are
        # always placed first.

        bb_verts = dict()
        bbs = sorted(
            mol.building_block_vertices,
            key=lambda bb: len(bb.func_groups),
            reverse=True
        )
        for bb in bbs:
            bb_verts[bb] = mol.building_block_vertices[bb]
        mol.building_block_vertices = bb_verts
        return super()._prepare(mol)

    def _clean_up(self, mol):
        mol.num_windows = self.num_windows
        mol.num_window_types = self.num_window_types
        return super()._clean_up(mol)

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

    def __repr__(self):
        vertex_alignments = ', '.join(
            f'{v.id}: {v.get_aligner_edge()}'
            for v in self.vertices
        )
        return (
            f'cage.{self.__class__.__name__}('
            f'vertex_alignments={{{vertex_alignments}}})'
        )
