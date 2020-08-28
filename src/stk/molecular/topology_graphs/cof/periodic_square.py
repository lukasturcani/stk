"""
Periodic Square
===============

"""

import numpy as np
from ...reactions import GenericReactionFactory
from .cof import Cof
from .vertices import _LinearCofVertex, _NonLinearCofVertex
from ..topology_graph import Edge
from ...periodic_info import PeriodicInfo


class PeriodicSquare(Cof):
    """
    Represents a periodic square COF topology graph.

    Building blocks with four and two functional groups are required
    for this topology graph.

    See :class:`.Cof` for more details and examples.

    """

    def __init__(
        self,
        building_blocks,
        lattice_size,
        vertex_alignments=None,
        reaction_factory=GenericReactionFactory(),
        num_processes=1,
    ):
        """
        Initialize a :class:`.PeriodicSquare` instance.

        Parameters
        ----------
        building_blocks : :class:`tuple` or :class:`dict`
            Can be a :class:`tuple` of :class:`.BuildingBlock`
            instances, which should be placed on the topology graph.

            Can also be a :class:`dict` which maps the
            :class:`.BuildingBlock` instances to the ids of the
            vertices it should be placed on. A :class:`dict` is
            required when there are multiple building blocks with the
            same number of functional groups, because in this case
            the desired placement is ambiguous.

        lattice_size : :class:`tuple` of :class:`int`
            The size of the lattice in the x, y and z directions.

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

        Raises
        ------
        :class:`AssertionError`
            If the any building block does not have a
            valid number of functional groups.

        :class:`ValueError`
            If the there are multiple building blocks with the
            same number of functional_groups in `building_blocks`,
            and they are not explicitly assigned to vertices. The
            desired placement of building blocks is ambiguous in
            this case.

        :class:`~.cof.UnoccupiedVertexError`
            If a vertex of the COF topology graph does not have a
            building block placed on it.

        :class:`~.cof.OverlyOccupiedVertexError`
            If a vertex of the COF topology graph has more than one
            building block placed on it.

        """

        super().__init__(
            building_blocks=building_blocks,
            lattice_size=lattice_size,
            periodic=True,
            vertex_alignments=vertex_alignments,
            reaction_factory=reaction_factory,
            num_processes=num_processes,
        )

    def get_periodic_info(self):
        """
        Get unit cell information of periodic topology graph.

        Returns
        -------
        :class:`.PeriodicInfo`
            Periodic cell information.

        """

        lattice_constants = self._get_lattice_constants()

        return PeriodicInfo(
            vector_1=(
                lattice_constants[0]*self._lattice_size[0]*self._scale
            ),
            vector_2=(
                lattice_constants[1]*self._lattice_size[1]*self._scale
            ),
            vector_3=(
                lattice_constants[2]*self._lattice_size[2]*self._scale
            ),
        )

    _lattice_constants = _a, _b, _c = (
        np.array([1., 0., 0.]),
        np.array([0., 1., 0.]),
        np.array([0., 0., 1.])
    )

    _vertex_prototypes = (
        _NonLinearCofVertex(0, (0.5)*_a + (0.5)*_b + (0.5)*_c),
    )
    _vertex_prototypes = (
        *_vertex_prototypes,
        _LinearCofVertex.init_at_shifted_center(
            id=1,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[0]),
            cell_shifts=((0, 0, 0), (1, 0, 0)),
            lattice_constants=_lattice_constants,
        ),
        _LinearCofVertex.init_at_shifted_center(
            id=2,
            vertices=(_vertex_prototypes[0], _vertex_prototypes[0]),
            cell_shifts=((0, 0, 0), (0, 1, 0)),
            lattice_constants=_lattice_constants,
        ),

    )

    _edge_prototypes = (
        Edge(0, _vertex_prototypes[1], _vertex_prototypes[0]),
        Edge(
            id=1,
            vertex1=_vertex_prototypes[1],
            vertex2=_vertex_prototypes[0],
            periodicity=(1, 0, 0),
        ),
        Edge(2, _vertex_prototypes[2], _vertex_prototypes[0]),
        Edge(
            id=3,
            vertex1=_vertex_prototypes[2],
            vertex2=_vertex_prototypes[0],
            periodicity=(0, 1, 0),
        ),
    )
