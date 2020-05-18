"""
M4L4 Square
===========

"""

from ..cage import Cage
from ..vertices import _LinearCageVertex
from ...topology_graph import Edge
from ....reactions import GenericReactionFactory


class M4L4Square(Cage):
    """
    Represents a cage topology graph.

    Both `corner` and `linker` vertices require building blocks with
    two functional groups for this topology.

    When using a :class:`dict` for initialization, a
    :class:`.BuildingBlock` needs to be assigned to each of the
    following numbers:

        | corners: (0, 1, 2, 3)
        | linkers: (4, 5, 6, 7)

    See :class:`.Cage` for more details and examples.

    """

    def __init__(
        self,
        corners,
        linkers,
        vertex_alignments=None,
        reaction_factory=GenericReactionFactory(),
        num_processes=1,
    ):
        """
        Initialize a :class:`.M4L4Square`.

        Parameters
        ----------
        corners : :class:`dict` or :class:`.BuildingBlock`
            Can be a :class:`dict` which maps the
            :class:`.BuildingBlock` instances to the ids of the
            vertices it should be placed on.

            Can also be a :class:`.BuildingBlock` instance, which
            should be placed on all corner vertices on the topology
            graph.

        linkers : :class:`dict` or :class:`.BuildingBlock`
            Can be a :class:`dict` which maps the
            :class:`.BuildingBlock` instances to the ids of the
            vertices it should be placed on.

            Can also be a :class:`.BuildingBlock` instance, which
            should be placed on all linker vertices on the topology
            graph.

        vertex_alignments : :class:`dict`, optional
            A mapping from the id of a :class:`.Vertex`
            to an :class:`.Edge` connected to it.
            The :class:`.Edge` is used to align the first
            :class:`.FunctionalGroup` of a :class:`.BuildingBlock`
            placed on that vertex. Only vertices which need to have
            their default edge changed need to be present in the
            :class:`dict`. If ``None`` then the default edge is used
            for each vertex. Changing which :class:`.Edge` is used will
            mean that the topology graph represents different
            structural isomers. The edge is referred to by a number
            between ``0`` (inclusive) and the number of edges the
            vertex is connected to (exclusive).

        reaction_factory : :class:`.ReactionFactory`, optional
            The reaction factory to use for creating bonds between
            building blocks.

        num_processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        """

        if isinstance(corners, dict):
            building_blocks = corners
        else:
            building_blocks = {corners: (0, 1, 2, 3)}

        if isinstance(linkers, dict):
            linkers_dict = linkers
        else:
            linkers_dict = {linkers: (4, 5, 6, 7)}

        building_blocks.update(
            (building_block, vertices)
            for building_block, vertices in linkers_dict.items()
        )

        super().__init__(
            building_blocks,
            vertex_alignments=vertex_alignments,
            reaction_factory=reaction_factory,
            num_processes=num_processes,
        )

    _vertex_prototypes = (
        _LinearCageVertex(0, [1, 1, 0]),
        _LinearCageVertex(1, [1, -1, 0]),
        _LinearCageVertex(2, [-1, -1, 0]),
        _LinearCageVertex(3, [-1, 1, 0]),

        _LinearCageVertex(4, [1, 0, 0], False),
        _LinearCageVertex(5, [0, -1, 0], False),
        _LinearCageVertex(6, [-1, 0, 0], False),
        _LinearCageVertex(7, [0, 1, 0], False),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[4]),
        Edge(1, _vertex_prototypes[1], _vertex_prototypes[4]),

        Edge(2, _vertex_prototypes[1], _vertex_prototypes[5]),
        Edge(3, _vertex_prototypes[2], _vertex_prototypes[5]),

        Edge(4, _vertex_prototypes[2], _vertex_prototypes[6]),
        Edge(5, _vertex_prototypes[3], _vertex_prototypes[6]),

        Edge(6, _vertex_prototypes[3], _vertex_prototypes[7]),
        Edge(7, _vertex_prototypes[0], _vertex_prototypes[7]),
    )

    _num_windows = 1
    _num_window_types = 1
