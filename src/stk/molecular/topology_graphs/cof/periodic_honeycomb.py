"""
Periodic Honeycomb
==================

"""

from collections import Counter
from .cof import OverlyOccupiedVertexError, UnoccupiedVertexError
from .honeycomb import Honeycomb
from ...reactions import GenericReactionFactory
from ..topology_graph import TopologyGraph


class PeriodicHoneycomb(TopologyGraph):
    """
    Represents a periodic honeycomb COF topology graph.

    Building blocks with three and two functional groups are required
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

        self._internal = Honeycomb(
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
        """
        Yield the building blocks.

        Building blocks are yielded in an order based on their
        position in the topology graph. For two equivalent
        topology graphs, but with different building blocks,
        equivalently positioned building blocks will be yielded at the
        same time.

        Yields
        ------
        :class:`.BuildingBlock`
            A building block of the topology graph.

        """

        for building_block in self._internal.get_building_blocks():
            yield building_block

    def get_num_building_block(self, building_block):
        """
        Get the number of times `building_block` is present.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block whose frequency in the topology graph
            is desired.

        Returns
        -------
        :class:`int`
            The number of times `building_block` is present in the
            topology graph.

        """

        return self._internal.get_num_building_block(building_block)

    @classmethod
    def _check_building_block_vertices(
        cls,
        num_vertices,
        building_block_vertices,
    ):
        unassigned_ids = set(range(num_vertices))
        assigned_ids = set()
        vertices = (
            vertex
            for vertices_ in building_block_vertices.values()
            for vertex in vertices_
        )
        for vertex in vertices:
            if vertex.get_id() in assigned_ids:
                raise OverlyOccupiedVertexError(
                    f'Vertex {vertex.get_id()} has multiple building '
                    'blocks placed on it.'
                )
            assigned_ids.add(vertex.get_id())
            unassigned_ids.remove(vertex.get_id())

        if unassigned_ids:
            raise UnoccupiedVertexError(
                'The following vertices are unoccupied '
                f'{unassigned_ids}.'
            )

    def clone(self):

        return self._internal.clone()

    def _get_edge_groups(self, edges):
        """
        Get the edge groups for the COF.

        Parameters
        ----------
        edges : :class:`iterable` of :class:`.Edge`
            The edges which are to be placed into edge groups.

        Returns
        -------
        :class:`tuple` of :class:`.EdgeGroup`
            The edge groups. These are returned when the lattice is
            not periodic, which means that some edges should not be
            part of an edge group, in order to prevent bond formation
            across the edges of the lattice.

        None : :class:`NoneType`
            If an edge group is to be created for every single edge.

        """

        return self._internal._get_edge_groups()

    def _get_vertices(self, lattice):
        """
        Get the vertices in the `lattice`.

        Parameters
        ----------
        lattice : :class:`list`
            A nested list which can be in the form
            ``lattice[x][y][z][vertex_id]`` which returns a vertex
            in the (x, y, z) cell of the lattice.

        Returns
        -------
        :class:`tuple` of :class:`.Vertex`
            All the vertices extracted from `lattice`.

        """

        return self._internal._get_vertices()

    def _get_lattice(self, vertex_alignments):
        """
        Get the vertices of the topology graph instance.

        Parameters
        ---------
        vertex_alignments : :class:`dict`
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

        Returns
        -------
        :class:`list`
            A nested :class:`list` which can be indexed as
            ``vertices[x][y][z]``, which will return a :class:`dict`
            for the unit cell at (x, y, z). The :class:`dict` maps
            the vertices in :attr:`_vertex_prototypes` to its clone for
            that unit cell.

        """

        return self._internal._get_lattice()

    def _get_edges(self, lattice):
        """
        Create the edges of the topology graph instance.

        Parameters
        ----------
        lattice : :class:`list`
            A nested :class:`list` which can be indexed as
            ``vertices[x][y][z]``, which will return a :class:`dict`
            for the unit cell at (x, y, z). The :class:`dict` maps
            the vertices in :attr:`_vertex_prototypes` to the clones
            for that unit cell.

        Returns
        -------
        :class:`tuple` of :class:`.Edge`
            The edges of the topology graph instance.

        """

        return self._internal._get_edges()

    @classmethod
    def _get_building_block_vertices(
        cls,
        building_blocks,
        vertices,
        edges,
    ):
        """
        Map building blocks to the vertices of the graph.

        Parameters
        ----------
        building_blocks : :class:`iterable` of :class:`.BuildingBlock`
            The building blocks which need to be mapped to `vertices`.

        vertices : :class:`iterable` of :class:`.Vertex`
            The vertices which need to have a building block map to
            them.

        edges : :class:`iterable`of :class:`.Edge`
            The edges of the graph.

        Returns
        -------
        :class:`dict`
            Maps each building block in `building_blocks` to a
            :class:`list` of the :class:`.Vertex` instances it should
            be placed on.

        Raises
        ------
        :class:`AssertionError`
            If the any building block does not have a
            valid number of functional groups.

        :class:`ValueError`
            If there are multiple building blocks with the same number
            of functional groups.

        """

        building_blocks_by_degree = {}
        for building_block in building_blocks:
            num_fgs = building_block.get_num_functional_groups()
            assert (
                num_fgs in cls._allowed_degrees
            ), (
                'The number of functional groups in '
                f'{building_block} needs to be one of '
                f'{tuple(cls._allowed_degrees)}, but is '
                'currently '
                f'{building_block.get_num_functional_groups()}.'
            )
            if num_fgs in building_blocks_by_degree:
                raise ValueError(
                    'If there are multiple building blocks with the '
                    'same number of functional groups, '
                    'building_block_vertices must be set explicitly.'
                )
            building_blocks_by_degree[num_fgs] = building_block

        vertex_degrees = Counter(
            vertex_id
            for edge in edges
            for vertex_id in edge.get_vertex_ids()
        )

        building_block_vertices = {}
        for vertex in vertices:
            vertex_degree = vertex_degrees[vertex.get_id()]
            building_block = building_blocks_by_degree[vertex_degree]
            building_block_vertices[building_block] = (
                building_block_vertices.get(building_block, [])
            )
            building_block_vertices[building_block].append(vertex)
        return building_block_vertices

    def _get_lattice_constants(self):
        return self._internal._get_lattice_constants()

    def _get_scale(self, building_block_vertices):
        return self._internal._get_scale(building_block_vertices)

    def get_periodic_cell(self):
        """
        Get unit cell matrix of periodic topology graph.

        Returns
        -------
        cell_matrix : :class:`tuple` of :class:`np.array`
            Tuple of cell lattice vectors (shape: (3,)) in Angstrom.

        """

        lattice_constants = self._get_lattice_constants()
        cell_matrix = tuple(
            i*j*self._internal._scale
            for i, j in zip(
                lattice_constants,
                self._internal._lattice_size
            )
        )

        return cell_matrix

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
