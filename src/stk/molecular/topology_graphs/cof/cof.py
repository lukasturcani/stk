import itertools as it
from collections import Counter
import numpy as np
from functools import partial
from operator import getitem

from ..topology_graph import TopologyGraph, EdgeGroup
from .edge import _CofEdge
from ...reactions import GenericReactionFactory


class Cof(TopologyGraph):
    def __init__(
        self,
        building_blocks,
        lattice_size,
        periodic=False,
        vertex_alignments=None,
        reaction_factory=GenericReactionFactory(),
        num_processes=1,
    ):

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
            get_vertex = partial(getitem, vertices)
            building_block_vertices = {
                building_block: map(get_vertex, ids)
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

        super().__init__(
            building_block_vertices=building_block_vertices,
            edges=edges,
            reaction_factory=reaction_factory,
            construction_stages=(),
            num_processes=num_processes,
            edge_groups=self._get_edge_groups(edges),
        )

    def _get_edge_groups(self, edges):
        if self._periodic:
            return None

        return tuple(
            EdgeGroup((edge,)) for edge in edges
            if not edge.is_periodic()
        )

    def _get_vertices(self, lattice):
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
                    {} for k in zdim
                ]
                for j in ydim
            ]
            for i in xdim
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

    @staticmethod
    def _get_building_block_vertices(building_blocks, vertices, edges):
        bb_by_degree = {}
        for bb in building_blocks:
            if bb.get_num_functional_groups() in bb_by_degree:
                raise ValueError(
                    'If there are multiple building blocks with the '
                    'same number of functional groups, '
                    'building_block_vertices must be set explicitly.'
                )
            bb_by_degree[bb.get_num_functional_groups()] = bb

        vertex_degrees = Counter(
            vertex_id
            for edge in edges
            for vertex_id in (
                edge.get_vertex1_id(),
                edge.get_vertex2_id(),
            )
        )

        building_block_vertices = {}
        for vertex in vertices:
            bb = bb_by_degree[vertex_degrees[vertex.get_id()]]
            building_block_vertices[bb] = (
                building_block_vertices.get(bb, [])
            )
            building_block_vertices[bb].append(vertex)
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
