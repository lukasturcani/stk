"""
Covalent Organic Framework
==========================

.. toctree::
    :maxdepth: 2

    Hexagonal <stk.molecular.topology_graphs.cof.hexagonal>
    Honeycomb <stk.molecular.topology_graphs.cof.honeycomb>
    Kagome <stk.molecular.topology_graphs.cof.kagome>
    Linkerless Honeycomb <\
stk.molecular.topology_graphs.cof.linkerless_honeycomb\
>
    Square <stk.molecular.topology_graphs.cof.square>

"""

import itertools as it
from collections import Counter
import numpy as np
from functools import partial
from operator import getitem

from ..topology_graph import TopologyGraph, EdgeGroup
from .edge import _CofEdge
from ...reactions import GenericReactionFactory


class UnoccupiedVertexError(Exception):
    """
    When a COF vertex is not occupied by a building block.

    """

    pass


class OverlyOccupiedVertexError(Exception):
    """
    When a COF vertex is occupied by more than one building block.

    """

    pass


class Cof(TopologyGraph):
    """
    Represents a COF topology graph.

    Notes
    -----
    COF topologies are added by creating a subclass, which defines
    the :attr:`_vertex_prototypes` and :attr:`_edge_prototypes` class
    attributes.

    Examples
    --------
    *Subclass Implementation*

    The source code of the subclasses, listed in :mod:`~.cof.cof`,
    can serve as good examples.

    *Basic Construction*

    :class:`.Cof` instances can be made by providing the building
    block molecules and lattice size only (using :class:`.Honeycomb`
    as an example)

    .. code-block:: python

        import stk

        bb1 = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])
        bb2 = stk.BuildingBlock('BrCC(CBr)CBr', [stk.BromoFactory()])
        cof = stk.ConstructedMolecule(
            topology_graph=stk.cof.Honeycomb((bb1, bb2), (3, 3, 1)),
        )

    *Structural Isomer Construction*

    Different structural isomers of COFs can be made by using the
    `vertex_alignments` optional parameter

    .. code-block:: python

        cof2 = stk.ConstructedMolecule(
            topology_graph=stk.cof.Honeycomb(
                building_blocks=(bb1, bb2),
                lattice_size=(3, 3, 1),
                vertex_alignments={0: 2, 1: 1, 2: 1},
            ),
        )

    The parameter maps the id of a vertex to a number
    between 0 (inclusive) and the number of edges the vertex is
    connected to (exclusive). So a vertex connected to three edges
    can be mapped to ``0``, ``1`` or ``2``.

    By changing which edge each vertex is aligned with, a different
    structural isomer of the COF can be formed.

    *Multi-Building Block COF Construction*

    You can also build COFs with multiple building blocks, but,
    if you have multiple building blocks with the same number of
    functional groups, you have to assign each building block to the
    vertex you want to place it on

    .. code-block:: python

        bb1 = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])
        bb2 = stk.BuildingBlock('BrCNCBr', [stk.BromoFactory()])
        bb3 = stk.BuildingBlock('BrCC(CBr)CBr', [stk.BromoFactory()])
        bb4 = stk.BuildingBlock('BrCC(NCBr)CBr', [stk.BromoFactory()])

        cof = stk.ConstructedMolecule(
            topology_graph=stk.cof.Honeycomb(
                # building_blocks is now a dict, which maps building
                # blocks to the id of the vertices it should be placed
                # on. You can use ranges to specify the ids.
                building_blocks={
                    bb1: (2, 4, 7, 8, 9, 12, 13, 14, 17, 18, 19),
                    bb2: 3,
                    bb3: (0, 10, 11, 15, 16),
                    bb4: (1, 5, 6),
                },
                lattice_size=(2, 2, 1),
            ),
        )

    You can combine this with the `vertex_alignments` parameter

    .. code-block:: python

        cof2 = stk.ConstructedMolecule(
            topology_graph=stk.cof.Honeycomb(
                building_blocks={
                    bb1: (2, 4, 7, 8, 9, 12, 13, 14, 17, 18, 19),
                    bb2: 3,
                    bb3: (0, 10, 11, 15, 16),
                    bb4: (1, 5, 6),
                },
                lattice_size=(2, 2, 1),
                vertex_alignments={0: 2, 1: 1, 2: 1},
            ),
        )

    """

    def __init_subclass__(cls, **kwargs):
        vertex_degrees = Counter(
            vertex_id
            for edge in cls._edge_prototypes
            for vertex_id in edge.get_vertex_ids()
        )
        cls._allowed_degrees = set(vertex_degrees.values())

    def __init__(
        self,
        building_blocks,
        lattice_size,
        periodic=False,
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

        periodic : :class:`bool`, optional
            Toggle the construction of a periodic molecule. If
            ``True``, periodic bonds will be made across the edges of
            the lattice.

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

        self._vertex_alignments = (
            dict(vertex_alignments)
            if vertex_alignments is not None
            else {}
        )
        self._lattice_size = lattice_size
        self._periodic = periodic

        lattice = self._get_lattice(self._vertex_alignments)
        edges = self._get_edges(lattice)
        vertices = self._get_vertices(lattice)

        if isinstance(building_blocks, dict):
            for building_block in building_blocks:
                assert (
                    building_block.get_num_functional_groups()
                    in self._allowed_degrees
                ), (
                    'The number of functional groups in '
                    f'{building_block} needs to be one of '
                    f'{tuple(self._allowed_degrees)}, but is '
                    'currently '
                    f'{building_block.get_num_functional_groups()}.'
                )
            get_vertex = partial(getitem, vertices)
            building_block_vertices = {
                building_block: tuple(map(
                    get_vertex,
                    # Account for the fact that a building block can
                    # be mapped to a single int.
                    (ids, ) if isinstance(ids, int) else ids
                ))
                for building_block, ids in building_blocks.items()
            }
        else:
            building_block_vertices = (
                self._get_building_block_vertices(
                    building_blocks=building_blocks,
                    vertices=vertices,
                    edges=edges,
                )
            )

        self._check_building_block_vertices(
            num_vertices=(
                np.product(lattice_size)*len(self._vertex_prototypes)
            ),
            building_block_vertices=building_block_vertices,
        )
        super().__init__(
            building_block_vertices=building_block_vertices,
            edges=edges,
            reaction_factory=reaction_factory,
            construction_stages=(),
            num_processes=num_processes,
            edge_groups=self._get_edge_groups(edges),
        )

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
        clone = super().clone()
        clone._vertex_alignments = dict(self._vertex_alignments)
        clone._lattice_size = self._lattice_size
        clone._periodic = self._periodic
        return clone

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

        if self._periodic:
            return None

        return tuple(
            EdgeGroup((edge,)) for edge in edges
            if not edge.is_periodic()
        )

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

        xdim, ydim, zdim = self._lattice_size
        num_vertices = xdim*ydim*zdim*len(self._vertex_prototypes)
        vertices = [None for i in range(num_vertices)]
        for x, y, z in it.product(
            range(xdim),
            range(ydim),
            range(zdim),
        ):
            for vertex in lattice[x][y][z].values():
                vertices[vertex.get_id()] = vertex
        return tuple(vertices)

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

        xdim, ydim, zdim = (range(dim) for dim in self._lattice_size)
        # vertex_clones is indexed as vertex_clones[x][y][z]
        lattice = [
            [
                [
                    {} for _ in zdim
                ]
                for _ in ydim
            ]
            for _ in xdim
        ]
        # Make a clone of each vertex for each unit cell.
        cells = it.product(xdim, ydim, zdim)
        vertices = it.product(cells, self._vertex_prototypes)
        for id_, (cell, vertex) in enumerate(vertices):
            x, y, z = cell
            shift = sum(
                axis * dim
                for axis, dim in zip(cell, self._lattice_constants)
            )
            lattice[x][y][z][vertex.get_id()] = vertex.__class__(
                id=id_,
                position=vertex.get_position() + shift,
                aligner_edge=vertex_alignments.get(id_, 0),
                cell=cell,
            )
        return lattice

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

        edge_clones = []
        # Make a clone for each edge for each unit cell.
        xdim, ydim, zdim = (range(dim) for dim in self._lattice_size)
        cells = it.product(xdim, ydim, zdim)
        edges = it.product(cells, self._edge_prototypes)
        for id_, (cell, edge) in enumerate(edges):
            x, y, z = cell
            # The cell in which the second vertex of the edge is found.
            periodic_cell = np.array(cell) + edge.get_periodicity()
            # Wrap around periodic cells, ie those that are less than 0
            # or greater than the lattice size along any dimension.
            dims = zip(periodic_cell, self._lattice_size)
            x2, y2, z2 = np.array([
                (dim+max_dim) % max_dim
                for dim, max_dim in dims
            ])
            # The edge is not periodic if periodic_cell did not
            # have to wrap around.
            dims = zip(periodic_cell, self._lattice_size)
            edge_is_not_periodic = all(
                dim >= 0 and dim < max_dim
                for dim, max_dim in dims
            )
            edge_clones.append(_CofEdge(
                parent_id=edge.get_id(),
                id=id_,
                vertex1=lattice[x][y][z][edge.get_vertex1_id()],
                vertex2=lattice[x2][y2][z2][edge.get_vertex2_id()],
                periodicity=(
                    (0, 0, 0)
                    if edge_is_not_periodic
                    else edge.get_periodicity()
                ),
            ))

        return tuple(edge_clones)

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
        return self._lattice_constants

    def _get_scale(self, building_block_vertices):
        return 5*max(
            bb.get_maximum_diameter()
            for bb in building_block_vertices
        )

    def __repr__(self):
        x, y, z = self._lattice_size
        periodic = ', periodic=True' if self._periodic else ''
        vertex_alignments = (
            f', vertex_alignments={self._vertex_alignments}'
            if self._vertex_alignments
            else ''
        )
        return (
            f'cof.{self.__class__.__name__}('
            f'({x}, {y}, {z})'
            f'{vertex_alignments}'
            f'{periodic})'
        )
