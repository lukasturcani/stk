import numpy as np
from scipy.spatial.distance import euclidean

from stk.utilities import (
    get_acute_vector,
    get_plane_normal,
    vector_angle,
    normalize_vector,
)
from ..topology_graph import Vertex


class _CageVertex(Vertex):
    """
    Represents a vertex of a :class:`.Cage`.

    """

    def __init__(
        self,
        id,
        position,
        use_neighbor_placement=True,
        aligner_edge=0,
    ):
        """
        Initialize a :class:`._CageVertex`.

        Parameters
        ----------
        id : :class:`int`
            The id of the vertex.

        position : :class:`tuple` of :class:`float`
            The position of the vertex.

        use_neighbor_placement : :class:`bool`, optional
            If ``True``, the position of the vertex will be updated
            based on the neighboring functional groups.

        aligner_edge : :class:`int`, optional
            The edge which is used to align the :class:`.BuildingBlock`
            placed on the vertex. The first :class:`.FunctionalGroup`
            is rotated such that it lies exactly on this
            :class:`.Edge`. Must be between ``0`` and the number of
            edges the vertex is connected to.

        """

        self._use_neighbor_placement = use_neighbor_placement
        self._aligner_edge = aligner_edge
        super().__init__(id, position)

    def clone(self):
        clone = super().clone()
        clone._aligner_edge = self._aligner_edge
        clone._use_neighbor_placement = self._use_neighbor_placement
        return clone

    def _with_aligner_edge(self, aligner_edge):
        """
        Modify the instance.

        """

        self._aligner_edge = aligner_edge
        return self

    def with_aligner_edge(self, aligner_edge):
        """
        Return a clone with a different `aligner_edge`.

        Parameters
        ----------
        aligner_edge : :class:`int`
            The aligner edge of the clone.

        Returns
        -------
        :class:`._CageVertex`
            The clone. Has the same type as the original instance.

        """

        return self.clone()._with_aligner_edge(aligner_edge)

    def use_neighbor_placement(self):
        """
        ``True`` if the position should be updated based on neighbors.

        Returns
        -------
        :class:`bool`
            ``True`` if the position of the vertex should be updated
            based on the positions of functional groups on neighboring
            vertices.

        """

        return self._use_neighbor_placement

    @classmethod
    def init_at_center(cls, id, vertices):
        """
        Initialize a :class:`._CageVertex` in the middle of `vertices`.

        Parameters
        ----------
        id : :class:`int`
            The id of the initialized vertex.

        vertices : :class:`tuple` of :class:`.Vertex`
            The vertices at whose center this one needs to be.

        Returns
        -------
        :class:`._CageVertex`
            The new vertex.

        """

        return cls(
            id=id,
            position=(
                sum(vertex.get_position() for vertex in vertices)
                / len(vertices)
            ),
        )

    def __str__(self):
        return (
            f'Vertex(id={self._id}, '
            f'position={self._position.tolist()}, '
            f'aligner_edge={self._aligner_edge})'
        )


class _LinearCageVertex(_CageVertex):
    def place_building_block(self, building_block, edges):
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
        fg_centroid = building_block.get_centroid(
            atom_ids=next(
                building_block.get_functional_groups()
            ).get_placer_ids(),
        )
        edge_position = edges[self._aligner_edge].get_position()
        edge_centroid = (
            sum(edge.get_position() for edge in edges) / len(edges)
        )
        building_block = building_block.with_rotation_between_vectors(
            start=fg_centroid - self._position,
            target=edge_position - edge_centroid,
            origin=self._position,
        )
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        return building_block.with_rotation_to_minimize_angle(
            start=core_centroid - self._position,
            target=self._position,
            axis=normalize_vector(
                edges[0].get_position() - edges[1].get_position()
            ),
            origin=self._position,
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        fg, = building_block.get_functional_groups(0)
        fg_position = building_block.get_centroid(fg.get_placer_ids())

        def fg_distance(edge):
            return euclidean(edge.get_position(), fg_position)

        edges = sorted(edges, key=fg_distance)
        return {
            fg_id: edge.get_id() for fg_id, edge in enumerate(edges)
        }


class _NonLinearCageVertex(_CageVertex):
    def place_building_block(self, building_block, edges):
        building_block = building_block.with_centroid(
            position=self._position,
            atom_ids=building_block.get_placer_ids(),
        )
        edge_centroid = (
            sum(edge.get_position() for edge in edges) / len(edges)
        )
        edge_normal = get_acute_vector(
            reference=edge_centroid,
            vector=get_plane_normal(
                points=np.array([
                    edge.get_position() for edge in edges
                ]),
            ),
        )
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        placer_centroid = building_block.get_centroid(
            atom_ids=building_block.get_placer_ids(),
        )
        building_block = building_block.with_rotation_between_vectors(
            start=get_acute_vector(
                reference=core_centroid - placer_centroid,
                vector=building_block.get_plane_normal(
                    atom_ids=building_block.get_placer_ids(),
                ),
            ),
            target=edge_normal,
            origin=self._position,
        )
        fg_bonder_centroid = building_block.get_centroid(
            atom_ids=next(
                building_block.get_functional_groups()
            ).get_placer_ids(),
        )
        edge_position = edges[self._aligner_edge].get_position()
        return building_block.with_rotation_to_minimize_angle(
            start=fg_bonder_centroid - self._position,
            target=edge_position - edge_centroid,
            axis=edge_normal,
            origin=self._position,
        ).get_position_matrix()

    def map_functional_groups_to_edges(self, building_block, edges):
        # The idea is to order the functional groups in building_block
        # by their angle with the vector running from the placer
        # centroid to fg0, going in the clockwise direction.
        # The edges are also ordered by their angle with the vector
        # running from the edge centroid to the aligner_edge,
        # going in the clockwise direction.
        #
        # Once the fgs and edges are ordered, zip and assign them.

        fg_sorter = _FunctionalGroupSorter(building_block)
        edge_sorter = _EdgeSorter(edges, fg_sorter.get_axis())
        return {
            fg_id: edge.get_id()
            for fg_id, edge in zip(
                fg_sorter.get_items(),
                edge_sorter.get_items(),
            )
        }


class _Sorter:
    """
    Sorts items according to their angle from a reference vector.

    """

    __slots__ = ['_items', '_reference', '_axis']

    def __init__(self, items, reference, axis):
        """
        Initialize a :class:`._Sorter`.

        Parameters
        ----------
        items : :class:`iterable` of :class:`object`
            The items to be sorted.

        reference : :class:`numpy.ndarray`
            The reference from which the angle is calculated.

        axis : :class:`numpy.ndarray`
            A vector orthogonal to `reference`, used to determine
            which direction is clockwise. Must be an immutable
            array.

        """

        self._items = items
        self._reference = reference
        self._axis = axis

    def _get_vector(self, item):
        """
        Get the vector according to which `item` should be sorted.

        Parameters
        ----------
        item : :class:`object`
            The item being sorted.

        Returns
        -------
        :class:`numpy.ndarray`
            The vector of the `item`, which should be used to measure
            its angle with respect to the *reference*.

        """

        raise NotImplementedError()

    def _get_angle(self, item):
        """
        Get the angle of `vector` relative to `reference`.

        Parameters
        ----------
        item : :class:`object`
            The item being sorted.

        Returns
        -------
        :class:`float`
            The angle between `item` and the reference vector.

        """

        vector = self._get_vector(item)
        theta = vector_angle(self._reference, vector)
        projection = vector @ self._axis
        if theta > 0 and projection < 0:
            return 2*np.pi - theta
        return theta

    def get_items(self):
        """
        Yield the sorted items.

        Yields
        ------
        :class:`object`
            An item.

        """

        yield from sorted(self._items, key=self._get_angle)

    def get_axis(self):
        """
        Get the axis used to determine which direction is clockwise.

        Returns
        -------
        :class:`numpy.ndarray`
            The axis. The array is immutable.

        """

        return self._axis


class _FunctionalGroupSorter(_Sorter):
    """
    Sorts functional groups according to their angle.

    """

    __slots__ = [
        '_items',
        '_reference',
        '_axis',
        '_placer_centroid',
    ]

    def __init__(self, building_block):
        """
        Initialize a :class:`._FunctionalGroupSorter` instance.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`

        """

        fg0_position = building_block.get_centroid(
            atom_ids=next(
                building_block.get_functional_groups()
            ).get_placer_ids(),
        )
        self._placer_centroid = placer_centroid = (
            building_block.get_centroid(
                atom_ids=building_block.get_placer_ids(),
            )
        )
        fg0_direction = fg0_position - placer_centroid
        core_centroid = building_block.get_centroid(
            atom_ids=building_block.get_core_atom_ids(),
        )
        axis = np.cross(
            fg0_direction,
            get_acute_vector(
                reference=core_centroid - placer_centroid,
                vector=building_block.get_plane_normal(),
            ),
        )
        axis.setflags(write=False)
        super().__init__(
            items=range(building_block.get_num_functional_groups()),
            reference=fg0_direction,
            axis=axis,
        )

    def _get_vector(self, item):
        building_block = self._building_block
        fg, = building_block.get_functional_groups(item)
        fg_position = building_block.get_centroid(fg.get_placer_ids())
        return fg_position - self._placer_centroid


class _EdgeSorter(_Sorter):
    """
    Sorted edges according to their angle.

    """

    __slots__ = [
        '_items',
        '_reference',
        '_axis',
        '_edge_centroid',
    ]

    def __init__(self, edges, aligner_edge, axis):
        """
        Initialize an :class:`._EdgeSorter` instance.

        Parameters
        ----------
        edges : :class:`iterable` of :class:`.Edge`
            The edges to sort.

        aligner_edge : :class:`.Edge`
            The edge in edges, used to calculate the reference vector.

        axis : :class:`numpy.ndarray`
            Must be immutable. The axis used to determine the clockwise
            direction.

        """

        self._edge_centroid = edge_centroid = (
            sum(edge.get_position() for edge in edges) / len(edges)
        )
        super().__init__(
            items=edges,
            axis=axis,
            reference=aligner_edge.get_position() - edge_centroid,
        )

    def _get_vector(self, item):
        return item.get_position() - self._edge_centroid
