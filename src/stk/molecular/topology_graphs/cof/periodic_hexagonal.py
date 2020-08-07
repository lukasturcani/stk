"""
Periodic Hexagonal
==================

"""

from .hexagonal import Hexagonal
from ...reactions import GenericReactionFactory
from ..topology_graph import TopologyGraph
from ...molecules import PeriodicInfo


class PeriodicHexagonal(TopologyGraph):
    """
    Represents a periodic hexagonal COF topology graph.

    Building blocks with six and two functional groups are required
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
        Initialize a :class:`.Cof` instance.

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

        self._internal = Hexagonal(
            building_blocks=building_blocks,
            lattice_size=lattice_size,
            periodic=True,
            vertex_alignments=vertex_alignments,
            reaction_factory=reaction_factory,
            num_processes=num_processes,
        )

    def construct(self):
        """
        Construct a :class:`.ConstructedMolecule`.

        Returns
        -------
        :class:`.ConstructionResult`
            The data describing the :class:`.ConstructedMolecule`.

        """

        return self._internal.construct()

    def get_building_blocks(self):
        yield from self._internal.get_building_blocks()

    def get_num_building_block(self, building_block):
        return self._internal.get_num_building_block(building_block)

    def clone(self):
        return self._internal.clone()

    def get_periodic_info(self):
        """
        Get unit cell information of periodic topology graph.

        Returns
        -------
        :class:`.PeriodicInfo`
            Class containing periodic information of cell.

        """

        lattice_constants = self._internal._get_lattice_constants()

        return PeriodicInfo(
            cell_matrix=tuple(
                i*j*self._internal._scale
                for i, j in zip(
                    lattice_constants,
                    self._internal._lattice_size
                )
            )
        )

    def __repr__(self):
        x, y, z = self._internal._lattice_size
        periodic = (
            ', periodic=True' if self._internal._periodic else ''
        )
        vertex_alignments = (
            f', vertex_alignments={self._internal._vertex_alignments}'
            if self._internal._vertex_alignments
            else ''
        )
        return (
            f'cof.{self.__class__.__name__}('
            f'({x}, {y}, {z})'
            f'{vertex_alignments}'
            f'{periodic})'
        )
