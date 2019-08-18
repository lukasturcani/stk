"""
Defines cages made from building blocks with 3 functional groups.

"""

import numpy as np

from .base import Cage,  _CageVertex
from ..topology_graph import Edge


class _OnePlusOneVertex(_CageVertex):
    def __init__(self, x, y, z, edge_normal, use_bonder_placement=True):
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

        edge_normal : :class:`list` of :class:`int`
            The edge plane normal to use.

        use_bonder_placement : :class:`bool`, optional
            If ``True``the position of the vertex will be updated such
            that it is in the middle of the neighboring bonder
            centroids, rather than in the middle of the neighboring
            vertices.

        """

        self._edge_normal = edge_normal
        super().__init__(x, y, z, use_bonder_placement)

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
        clone._edge_normal = list(self._edge_normal)
        return clone

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
        building_block.apply_rotation_between_vectors(
            start=building_block.get_bonder_plane_normal(),
            target=self._edge_normal,
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
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    _x = 1
    vertices = (
        _OnePlusOneVertex(_x, 0., 0., [1, 0, 0], False),
        _OnePlusOneVertex(-_x, 0., 0., [-1, 0, 0], False),

    )
    edges = (
        Edge(
            vertices[0], vertices[1], position=np.array([0., 1., 0.])
        ),
        Edge(
            vertices[0], vertices[1], position=np.array([0., -1., 1.])
        ),
        Edge(
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
        _CageVertex(_x, 0, -_x/np.sqrt(2), False),
        _CageVertex(-_x, 0, -_x/np.sqrt(2), False),
        _CageVertex(0, _x, _x/np.sqrt(2), False),
        _CageVertex(0, -_x, _x/np.sqrt(2), False)
    )

    edges = (
        Edge(vertices[0], vertices[1]),
        Edge(vertices[0], vertices[2]),
        Edge(vertices[0], vertices[3]),

        Edge(vertices[1], vertices[2]),
        Edge(vertices[1], vertices[3]),

        Edge(vertices[2], vertices[3])
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
        _CageVertex(-_x, _x, -_x, False),
        _CageVertex(-_x, -_x, -_x, False),
        _CageVertex(_x, _x, -_x, False),
        _CageVertex(_x, -_x, -_x, False),

        _CageVertex(-_x, _x, _x, False),
        _CageVertex(-_x, -_x, _x, False),
        _CageVertex(_x, _x, _x, False),
        _CageVertex(_x, -_x, _x, False)
    )

    edges = (
        Edge(vertices[0], vertices[1]),
        Edge(vertices[0], vertices[2]),
        Edge(vertices[0], vertices[4]),
        Edge(vertices[1], vertices[3]),
        Edge(vertices[1], vertices[5]),
        Edge(vertices[2], vertices[6]),
        Edge(vertices[2], vertices[3]),
        Edge(vertices[3], vertices[7]),
        Edge(vertices[4], vertices[6]),
        Edge(vertices[4], vertices[5]),
        Edge(vertices[5], vertices[7]),
        Edge(vertices[6], vertices[7])
    )

    num_windows = 6
    num_window_types = 1
