import numpy as np
from collections import defaultdict

from ..topology_graph import TopologyGraph, Vertex
from ....utilities import vector_theta


class _CageVertex(Vertex):
    """
    Represents a vertex of a :class:`.CageTopology`.

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

    def __init__(self, x, y, z):
        """
        Initialize a :class:`_CageVertex`.

        Parameters
        ----------
        x : :class:`float`
            The x coordinate.

        y : :class:`float`
            The y coordinate.

        z : :class:`float`
            The z coordinate.

        """

        # _neighbor_positions holds the bonder centroids of functional
        # groups on neighbor vertices connected to this vertex.
        self._neighbor_positions = []
        self.aligner_edge = None
        # id will be set automatically by CageTopology. This is because
        # _CageVertex is defined manually in a subclass of CageTopology
        # and writing the id for every vertex would be a pain.
        super().__init__(None, x, y, z)

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
        clone.aligner_edge = self.aligner_edge
        clone._neighbor_positions = list(self._neighbor_positions)
        return clone

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

        if len(self._neighbor_positions) == len(self.edges):
            self._update_position()

        if len(building_block.func_groups) == 2:
            return self._place_linear_building_block(building_block)
        return self._place_nonlinear_building_block(building_block)

    def _update_position(self):
        self._position = np.divide(
            np.sum(self._neighbor_positions, axis=0),
            len(self._neighbor_positions)
        )

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
        building_block.apply_rotation_to_minimize_theta(
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
        building_block.apply_rotation_to_minimize_theta(
            start=start,
            target=target,
            axis=edge_normal,
            origin=self._position
        )
        return building_block.get_position_matrix()

    def assign_func_groups_to_edges(self, building_block, fg_map):
        """
        Assign functional groups to edges.

        Each :class:`.FunctionalGroup` of the `building_block` needs
        to be associated with one of the :class:`.Edge` instances in
        :attr:`edges`. Then, using `fg_map`, the
        :class:`FunctionalGroup` instances in the molecule being
        constructed need to be assigned to those edges. This is
        because bonds need to be formed between functional groups of
        the molecule being constructed, not the `building_block`.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is needs to have
            functional groups assigned to edges.

        fg_map : :class:`dict`
            A mapping from :class:`.FunctionalGroup` instances in
            `building_block` to the equivalent
            :class:`.FunctionalGroup` instances in the molecule being
            constructed.

        Returns
        -------
        None : :class:`NoneType`

        """

        if len(building_block.func_groups) == 2:
            r = self._assign_func_groups_to_linear_edges(
                building_block=building_block,
                fg_map=fg_map
            )
        else:
            r = self._assign_func_groups_to_nonlinear_edges(
                building_block=building_block,
                fg_map=fg_map
            )

        bb_fgs = set(fg_map.values())

        for edge in self.edges:
            for vertex in edge.vertices:
                if vertex is self:
                    continue

                for func_group in edge.get_func_groups():
                    if func_group not in bb_fgs:
                        continue

                    vertex._neighbor_positions.append(
                        self._get_molecule_centroid(
                            atom_ids=func_group.get_bonder_ids()
                        )
                    )
        return r

    def _assign_func_groups_to_linear_edges(
        self,
        building_block,
        fg_map
    ):

        fg1, fg2 = sorted(
            building_block.func_groups,
            key=self._get_edge0_distance(building_block)
        )
        self.edges[0].assign_func_group(fg_map[fg1])
        self.edges[1].assign_func_group(fg_map[fg2])

    def _get_edge0_distance(self, building_block):
        aligner_coord = self.edges[0].get_position()

        def distance(fg):
            fg_coord = building_block.get_centroid(
                atom_ids=fg.get_bonder_ids()
            )
            displacement = aligner_coord - fg_coord
            return np.linalg.norm(displacement)

        return distance

    def _assign_func_groups_to_nonlinear_edges(
        self,
        building_block,
        fg_map
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
            building_block.func_groups,
            key=self._get_func_group_angle(
                building_block=building_block,
                fg0_direction=fg0_direction,
                bonder_centroid=bonder_centroid,
                axis=axis
            )
        )
        assert func_groups[0] is building_block.func_groups[0]

        edges = sorted(self.edges, key=self._get_edge_angle(axis))
        aligner_first = all(
            edges[0].get_position() == self.aligner_edge.get_position()
        )
        assert aligner_first

        for edge, func_group in zip(edges, func_groups):
            edge.assign_func_group(fg_map[func_group])

    @staticmethod
    def _get_func_group_angle(
        building_block,
        fg0_direction,
        bonder_centroid,
        axis
    ):

        def angle(func_group):
            coord = building_block.get_centroid(
                atom_ids=func_group.get_bonder_ids()
            )
            fg_direction = coord-bonder_centroid
            theta = vector_theta(fg0_direction, fg_direction)

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
            theta = vector_theta(
                vector1=edge_direction,
                vector2=aligner_edge_direction
            )

            projection = edge_direction @ axis
            if theta > 0 and projection < 0:
                return 2*np.pi - theta
            return theta

        return angle

    def __str__(self):
        x, y, z = self._position
        return (
            f'Vertex(id={self.id}, '
            f'position={[x, y, z]}, '
            f'aligner_edge={self.edges.index(self.aligner_edge)})'
        )


class CageTopology(TopologyGraph):
    """
    Represents a cage topology graph.

    Cage topologies are added by creating a subclass which defines the
    :attr:`vertices` and :attr:`edges` of the topology as class
    attributes.

    A :class:`CageTopology` subclass will add the attributes
    :attr:`num_windows` and :attr:`num_window_types` to each
    :class:`.ConstructedMolecule`.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    Examples
    --------
    :class:`CageTopology` instances can be made without supplying
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

        v0 = stk.FourPlusSix.vertices[0]
        v2 = stk.FourPlusSix.vertices[2]
        tetrahedron = stk.cage.FourPlusSix(
            vertex_alignments={
                v0: v0.edges[1],
                v2: v2.edges[2]
            }
        )
        cage2 = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topology_graph=tetrahedron
        )

    By changing which edge each vertex is aligned with, a different
    structural isomer of the cage can be formed.

    Note the in the `vertex_alignments` parameter the class vertices
    and edges are used, however when the `building_block_vertices`
    parameter is used, the instance vertices are used. **These are not
    interchangeable!**

    .. code-block:: python

        # Use the class vertices and edges to set vertex_alignments
        # and create a topology graph.
        v0 = stk.FourPlusSix.vertices[0]
        v2 = stk.FourPlusSix.vertices[2]
        tetrahedron = stk.cage.FourPlusSix(
            vertex_alignments={
                v0: v0.edges[1],
                v2: v2.edges[2]
            }
        )
        bb3 = stk.BuildingBlock('NCOCN', ['amine'])
        cage2 = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2, bb3],
            topology_graph=tetrahedron
            # Use the instance vertices in the building_block_vertices
            # parameter.
            building_block_vertices={
                bb1: tetrahedron.vertices[:2],
                bb2: tetrahedron.vertices[4:],
                bb3: tetrahedron.vertices[2:4]
            }
        )

    The example above also demonstrates how cages with many building
    blocks can be built. You can add as many :class:`.BuildingBlock`
    instances into `building_blocks` as you like. If you do not
    assign where each building block is placed with
    `building_block_vertices`, they will be placed on the
    :atttr:`vertices` of the :class:`.CageTopology` at random. Random
    placement will account for the fact that the length of
    :attr:`.BuildingBlock.func_groups` needs to match the number of
    edges connected to a vertex.

    """

    def __init_subclass__(cls, **kwargs):
        for i, vertex in enumerate(cls.vertices):
            vertex.id = i
        return super().__init_subclass__(**kwargs)

    def __init__(self, vertex_alignments=None, processes=1):
        """
        Initialize a :class:`.CageTopology`.

        Parmaeters
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

        processes : :class:`int`, optional
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
        vertex_alignments = _vertex_alignments

        vertex_clones = {}
        for vertex in self.vertices:
            vertex.aligner_edge = vertex_alignments.get(
                vertex,
                vertex.edges[0]
            )
            clone = vertex.clone(clear_edges=True)
            vertex_clones[vertex] = clone

        edge_clones = {}
        for edge in self.edges:
            edge_clones[edge] = edge.clone(vertex_clones)

        for vertex in vertex_clones.values():
            vertex.aligner_edge = edge_clones[vertex.aligner_edge]

        super().__init__(
            vertices=tuple(vertex_clones.values()),
            edges=tuple(edge_clones.values()),
            processes=processes
        )

    def _assign_building_blocks_to_vertices(
        self,
        mol,
        building_blocks
    ):
        """
        Assign `building_blocks` to :attr:`vertices`.

        This method will assign a random building block with the
        correct amount of functional groups to each vertex.

        Assignment is done by modifying
        :attr:`.ConstructedMolecule.building_block_vertices`.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The :class:`.ConstructedMolecule` instance being
            constructed.

        building_blocks : :class:`list` of :class:`.Molecule`
            The :class:`.BuildingBlock` and
            :class:`ConstructedMolecule` instances which
            represent the building block molecules used for
            construction. Only one instance is present per building
            block molecule, even if multiples of that building block
            join up to form the :class:`ConstructedMolecule`.

        Returns
        -------
        None : :class:`NoneType`

        """

        bb_by_degree = defaultdict(list)
        for bb in building_blocks:
            bb_by_degree[len(bb.func_groups)].append(bb)

        for vertex in self.vertices:
            bb = np.random.choice(bb_by_degree[len(vertex.edges)])
            mol.building_block_vertices[bb].append(vertex)

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
            f'{v.id}: {v.edges.index(v.aligner_edge)}'
            for v in self.vertices
        )

        return (
            f'cage.{self.__class__.__name__}('
            f'vertex_alignments={{{vertex_alignments}}})'
        )
