"""
Defines cages made from building blocks with 3 functional groups.

"""

import numpy as np

from .base import Cage,  _CageVertexData, _CageVertex
from ..topology_graph import EdgeData


class _OnePlusOneVertexData(_CageVertexData):
    def __init__(
        self,
        x,
        y,
        z,
        edge_normal,
        cell=None,
        aligner_edge=None,
        use_bonder_placement=True
    ):
        """
        Initialize a :class:`_OnePlusOneVertexData` instance.

        Parameters
        ----------
        x : :class:`float`
            The x coordinate.

        y : :class:`float`
            The y coordinate.

        z : :class:`float`
            The z coordinate.

        edge_normal : :class:`list` of :class:`int`
            The edge plane normal to use.

        cell : :class:`numpy.ndarray`, optional
            The unit cell in which the vertex is found.

        aligner_edge : :class:`int`, optional
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

        self.edge_normal = edge_normal
        super().__init__(
            x=x,
            y=y,
            z=z,
            cell=cell,
            aligner_edge=aligner_edge,
            use_bonder_placement=use_bonder_placement
        )

    def get_vertex(self):
        return _OnePlusOneVertex(self)


class _OnePlusOneVertex(_CageVertex):

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
        clone._edge_normal = list(self._edge_normal)
        return clone

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
        building_block.apply_rotation_between_vectors(
            start=building_block.get_bonder_plane_normal(),
            target=self._edge_normal,
            origin=self._position
        )
        fg_bonder_centroid = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_bonder_ids()
        )
        start = fg_bonder_centroid - self._position
        aligner_edge = edges[self._edge_ids[self.aligner_edge]]
        edge_coord = aligner_edge.get_position()
        connected_edges = tuple(edges[id_] for id_ in self._edge_ids)
        target = edge_coord - self._get_edge_centroid(connected_edges)
        building_block.apply_rotation_to_minimize_angle(
            start=start,
            target=target,
            axis=self._edge_normal,
            origin=self._position
        )
        return building_block.get_position_matrix()


class OnePlusOne(Cage):
    """
    Represents a capsule cage topology graph.

    See :class:`.Cage` for more details and examples.

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

    _x = 1
    vertices = (
        _OnePlusOneVertexData(_x, 0., 0., [1, 0, 0], False),
        _OnePlusOneVertexData(-_x, 0., 0., [-1, 0, 0], False),

    )
    edges = (
        EdgeData(
            vertices[0], vertices[1], position=np.array([0., 1., 0.])
        ),
        EdgeData(
            vertices[0], vertices[1], position=np.array([0., -1., 1.])
        ),
        EdgeData(
            vertices[0], vertices[1], position=np.array([0., -1., -1.])
        )
    )

    num_windows = 3
    num_window_types = 1


class TwoPlusTwo(Cage):
    """
    Represents a tetrahedron cage topology graph.

    See :class:`.Cage` for more details and examples.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    _x = 1
    vertices = (
        _CageVertexData(_x, 0, -_x/np.sqrt(2), False),
        _CageVertexData(-_x, 0, -_x/np.sqrt(2), False),
        _CageVertexData(0, _x, _x/np.sqrt(2), False),
        _CageVertexData(0, -_x, _x/np.sqrt(2), False)
    )

    edges = (
        EdgeData(vertices[0], vertices[1]),
        EdgeData(vertices[0], vertices[2]),
        EdgeData(vertices[0], vertices[3]),

        EdgeData(vertices[1], vertices[2]),
        EdgeData(vertices[1], vertices[3]),

        EdgeData(vertices[2], vertices[3])
    )

    num_windows = 4
    num_window_types = 1


class FourPlusFour(Cage):
    """
    Represents a cube cage topology graph.

    See :class:`.Cage` for more details and examples.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    _x = 1
    vertices = (
        _CageVertexData(-_x, _x, -_x, False),
        _CageVertexData(-_x, -_x, -_x, False),
        _CageVertexData(_x, _x, -_x, False),
        _CageVertexData(_x, -_x, -_x, False),

        _CageVertexData(-_x, _x, _x, False),
        _CageVertexData(-_x, -_x, _x, False),
        _CageVertexData(_x, _x, _x, False),
        _CageVertexData(_x, -_x, _x, False)
    )

    edges = (
        EdgeData(vertices[0], vertices[1]),
        EdgeData(vertices[0], vertices[2]),
        EdgeData(vertices[0], vertices[4]),
        EdgeData(vertices[1], vertices[3]),
        EdgeData(vertices[1], vertices[5]),
        EdgeData(vertices[2], vertices[6]),
        EdgeData(vertices[2], vertices[3]),
        EdgeData(vertices[3], vertices[7]),
        EdgeData(vertices[4], vertices[6]),
        EdgeData(vertices[4], vertices[5]),
        EdgeData(vertices[5], vertices[7]),
        EdgeData(vertices[6], vertices[7])
    )

    num_windows = 6
    num_window_types = 1
