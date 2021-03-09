"""
M3L3 Triangle
=============

"""

import numpy as np

from ..cage import Cage
from ..vertices import _LinearCageVertex
from ...topology_graph import Edge, NullOptimizer
from ....reactions import GenericReactionFactory


class M3L3Triangle(Cage):
    """
    Represents a cage topology graph.

    Both `corner` and `linker` vertices require building blocks with
    two functional groups for this topology. This class replaces the
    `building_blocks` parameter with the `corner` and `linker`
    parameters.

    See :class:`.Cage` for more details and examples.

    """

    def __init__(
        self,
        corners,
        linkers,
        vertex_alignments=None,
        reaction_factory=GenericReactionFactory(),
        num_processes=1,
        optimizer=NullOptimizer(),
    ):
        """
        Initialize a :class:`.M3L3Triangle`.

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

        optimizer : :class:`.Optimizer`, optional
            Used to optimize the structure of the constructed
            molecule.

        """

        if isinstance(corners, dict):
            building_blocks = corners
        else:
            building_blocks = {corners: (0, 1, 2)}

        if isinstance(linkers, dict):
            linkers_dict = linkers
        else:
            linkers_dict = {linkers: (3, 4, 5)}

        building_blocks.update(
            (building_block, vertices)
            for building_block, vertices in linkers_dict.items()
        )

        super().__init__(
            building_blocks,
            vertex_alignments=vertex_alignments,
            reaction_factory=reaction_factory,
            num_processes=num_processes,
            optimizer=optimizer,
        )

    _x = 2*np.sqrt(3)/4
    _y = 2
    _vertex_prototypes = (
        _LinearCageVertex(0, [0, _x, 0]),
        _LinearCageVertex(1, [_y/2, -_x, 0]),
        _LinearCageVertex(2, [-_y/2, -_x, 0]),
    )

    _vertex_prototypes = (
        *_vertex_prototypes,

        _LinearCageVertex.init_at_center(
            id=3,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[1]),
        ),
        _LinearCageVertex.init_at_center(
            id=4,
            vertices=(_vertex_prototypes[1], _vertex_prototypes[2]),
        ),
        _LinearCageVertex.init_at_center(
            id=5,
            vertices=(_vertex_prototypes[2], _vertex_prototypes[0]),
        ),
    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[0], _vertex_prototypes[3]),
        Edge(1, _vertex_prototypes[1], _vertex_prototypes[3]),

        Edge(2, _vertex_prototypes[1], _vertex_prototypes[4]),
        Edge(3, _vertex_prototypes[2], _vertex_prototypes[4]),

        Edge(4, _vertex_prototypes[2], _vertex_prototypes[5]),
        Edge(5, _vertex_prototypes[0], _vertex_prototypes[5]),
    )

    _num_windows = 1
    _num_window_types = 1
