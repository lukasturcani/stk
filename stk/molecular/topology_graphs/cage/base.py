import numpy as np
from collections import defaultdict

from ..topology_graph import TopologyGraph, Vertex, Edge
from ....utilities import vector_theta


class _CageVertex(Vertex):
    """
    Represents a vertex of a :class:`.CageTopology`.

    Attributes
    ----------
    edges : :class:`list` of :class:`.Edge`
        The edges the :class:`Vertex` is connected to.

    aligner_edge : :class:`.Edge`
        The :class:`.Edge` in :attr:`edges`, which is used to align the
        :class:`.BuildingBlock` placed on the vertex. The first
        :class:`.FunctionalGroup` in :attr:`.BuildingBlock.func_groups`
        is rotated such that it lies exactly on this :class:`.Edge`.

    """

    def __init__(self, x, y, z):
        self.aligner_edge = None
        super().__init__(x, y, z)

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

        if len(building_block.func_groups) == 2:
            return self._place_linear_building_block(building_block)
        return self._place_nonlinear_building_block(building_block)

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
            position=self._coord,
            atom_ids=building_block.get_bonder_ids()
        )
        fg_centroid = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_bonder_ids()
        )
        start = fg_centroid - self._coord
        edge_coord = self.aligner_edge.get_position()
        target = edge_coord - self._get_edge_centroid()
        building_block.apply_rotation_between_vectors(
            start=start,
            target=target,
            origin=self._coord
        )
        start = building_block.get_centroid_centroid_direction_vector()
        e0_coord = self.edges[0].get_position()
        e1_coord = self.edges[1].get_position()
        building_block.apply_rotation_to_minimize_theta(
            start=start,
            target=self._coord,
            axis=e0_coord-e1_coord,
            origin=self._coord,
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
            position=self._coord,
            atom_ids=building_block.get_bonder_ids()
        )
        edge_normal = self._get_edge_plane_normal(
            reference=self._get_edge_centroid()
        )
        building_block.apply_rotation_between_vectors(
            start=building_block.get_bonder_plane_normal(),
            target=edge_normal,
            origin=self._coord
        )
        fg_centroid = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_bonder_ids()
        )
        start = fg_centroid - self._coord
        edge_coord = self.aligner_edge.get_position()
        target = edge_coord - self._get_edge_centroid()
        building_block.apply_rotation_to_minimize_theta(
            start=start,
            target=target,
            axis=edge_normal,
            origin=self._coord
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
            return self._assign_func_groups_to_edges_linear(
                building_block=building_block,
                fg_map=fg_map
            )
        return self._assign_func_groups_to_edges_nonlinear(
            building_block=building_block,
            fg_map=fg_map
        )

    def _assign_func_groups_to_edges_linear(
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

    def _assign_func_groups_to_edges_nonlinear(
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
        func_groups = sorted(
            building_block.func_groups,
            key=self._get_func_group_angle(
                building_block=building_block,
                fg0_coord=fg0_coord,
                bonder_centroid=bonder_centroid
            )
        )
        edges = sorted(self.edges, key=self._get_edge_angle())

        for func_group, edge in zip(func_groups, edges):
            edge.assign_func_group(fg_map[func_group])

    @staticmethod
    def _get_func_group_angle(
        building_block,
        fg0_coord,
        bonder_centroid
    ):

        # This axis is used to figure out the clockwise direction.
        axis = np.cross(
            fg0_coord-bonder_centroid,
            building_block.get_bonder_plane_normal()
        )

        def angle(func_group):
            coord = building_block.get_centroid(
                atom_ids=func_group.get_bonder_ids()
            )
            theta = vector_theta(coord, fg0_coord)

            projection = coord @ axis
            if projection < 0:
                return 2*np.pi - theta
            return theta

        return angle

    def _get_edge_angle(self):

        aligner_edge_coord = self.aligner_edge.get_position()
        edge_centroid = self._get_edge_centroid()
        # This axis is used to figure out the clockwise direction.
        axis = np.cross(
            aligner_edge_coord-edge_centroid,
            self._get_edge_plane_normal(self._get_edge_centroid())
        )

        def angle(edge):
            coord = edge.get_position()
            theta = vector_theta(coord, aligner_edge_coord)

            projection = coord @ axis
            if projection < 0:
                return 2*np.pi - theta
            return theta

        return angle


class CageTopology(TopologyGraph):
    """
    Represents a cage topology graph.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        A class attribute. It holds vertices used to make a specific
        cage topology graph. This needs to be defined by a subclass.

    edges : :class:`tuple` of :class:`.Edge`
        A class attribute. It hold the edges used to make a specific
        cage topology graph. This needs to be defined by a subclass.

    """

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

        processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        """

        if vertex_alignments is None:
            vertex_alignments = {}

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
            vertices = [
                vertex_clones[vertex] for vertex in edge.vertices
            ]
            edge_clones[edge] = Edge(*vertices)

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
