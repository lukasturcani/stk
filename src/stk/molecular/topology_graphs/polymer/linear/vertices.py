import logging

from ...topology_graph import Vertex


logger = logging.getLogger(__name__)


class _LinearVertex(Vertex):
    """
    Represents a vertex in the middle of a linear polymer chain.

    """

    def __init__(self, id, position, flip):
        """
        Initialize a :class:`._LinearVertex` instance.

        Parameters
        ----------
        id : :class:`int`
            The id of the vertex.

        position : :class:`numpy.ndarray`
            The position of the vertex.

        flip : :class:`bool`
            If ``True`` any building block placed by the vertex will
            have its orientation along the chain flipped.

        """

        super().__init__(id, position)
        self._flip = flip

    def clone(self):
        clone = super().clone()
        clone._flip = self._flip
        return clone

    def place_building_block(self, building_block):
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
        fg1, fg2 = building_block.get_functional_groups()
        fg1_position = building_block.get_centroid(
            atom_ids=fg1.get_placer_ids(),
        )
        fg2_position = building_block.get_centroid(
            atom_ids=fg2.get_placer_ids(),
        )
        fg_displacement = fg1_position - fg2_position
        building_block = building_block.with_rotation_between_vectors(
            start=fg_displacement,
            target=[-1 if self._flip else 1, 0, 0],
            origin=self._position,
        )
        return building_block.get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        fg1_id, fg2_id = self._sort_functional_groups(building_block)
        edge1_id, edge2_id = self._sort_edges(edges)
        return {
            fg1_id: edge1_id,
            fg2_id: edge2_id,
        }

    @staticmethod
    def _sort_functional_groups(building_block):
        fg1, fg2 = building_block.get_functional_groups()
        x1, y1, z1 = building_block.get_centroid(
            atom_ids=fg1.get_placer_ids(),
        )
        x2, y2, z2 = building_block.get_centroid(
            atom_ids=fg2.get_placer_ids(),
        )
        return 0, 1 if x1 < x2 else 1, 0

    @staticmethod
    def _sort_edges(edges):
        edge1, edge2 = edges
        x1, y1, z1 = edge1.get_position()
        x2, y2, z2 = edge2.get_position()
        if x1 < x2:
            return edge1.get_id(), edge2.get_id()
        else:
            return edge2.get_id(), edge1.get_id()

    def __str__(self):
        return (
            f'Vertex(id={self.id}, '
            f'position={self._position.tolist()}, '
            f'flip={self._flip})'
        )


class _TerminalVertex(_LinearVertex):
    """
    Represents a vertex at the end of a polymer chain.

    Do not instantiate this class directly, use :class:`.HeadVertex` or
    :class:`.TailVertex` instead.

    Attributes
    ----------
    id : :class:`int`
        The id of the vertex. This should be its index in
        :attr:`TopologyGraph.vertices`.

    """

    def place_building_block(self, building_block, vertices, edges):
        """
        Place `building_block` on the :class:`.Vertex`.

        `building_block` is placed such that the centroid-centroid
        direction vector is always pointing away from the center of the
        chain, when `building_block` has only one
        :class:`.FunctionalGroup`.

        If the `building_block` has more than one
        :class:`.FunctionalGroup` then this method behaves in the same
        was as for :class:`.LinearVertex`.

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
            placed on the :class:`.Vertex`.

        """

        if len(building_block.func_groups) != 1:
            return super().place_building_block(
                building_block=building_block,
                vertices=vertices,
                edges=edges
            )

        building_block.set_centroid(
            position=self._position,
            atom_ids=building_block.get_bonder_ids(fg_ids=(0, )),
        )
        centroid_vector = (
            building_block.get_centroid_centroid_direction_vector()
        )
        building_block.apply_rotation_between_vectors(
            start=centroid_vector,
            # _cap_direction is defined by a subclass.
            target=[self._cap_direction, 0, 0],
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

        Raises
        ------
        :class:`ValueError`
            If `building_block` does not have one or two functional
            groups.

        """

        if len(building_block.func_groups) == 2:
            func_groups = building_block.func_groups
            fgs = sorted(
                range(len(func_groups)),
                key=lambda fg_id: building_block.get_centroid(
                    atom_ids=func_groups[fg_id].get_bonder_ids()
                )[0]
            )
            fg_index = 0 if self._cap_direction == 1 else -1
            return {fgs[fg_index]: self._edge_ids[0]}

        elif len(building_block.func_groups) == 1:
            return {0: self._edge_ids[0]}

        else:
            raise ValueError(
                'The building block of a polymer '
                'must have 1 or 2 functional groups.'
            )


class _HeadVertex(_TerminalVertex):
    """
    Represents a vertex at the head of a polymer chain.

    Attributes
    ----------
    id : :class:`int`
        The id of the vertex. This should be its index in
        :attr:`TopologyGraph.vertices`.

    """

    # The direction to use if the building block placed on the
    # vertex only has 1 FunctionalGroup.
    _cap_direction = -1


class _TailVertex(_TerminalVertex):
    """
    Represents a vertex at the tail of a polymer chain.

    Attributes
    ----------
    id : :class:`int`
        The id of the vertex. This should be its index in
        :attr:`TopologyGraph.vertices`.

    """

    # The direction to use if the building block placed on the
    # vertex only has 1 FunctionalGroup.
    _cap_direction = 1


