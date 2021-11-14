"""
Graph State
===========

"""

from collections import defaultdict

import numpy as np


class _GraphState:
    """
    The topology graph of a molecule under construction.

    """

    __slots__ = [
        '_vertex_building_blocks',
        '_vertices',
        '_edges',
        '_lattice_constants',
        '_vertex_edges',
        '_num_building_blocks',
    ]

    def __init__(
        self,
        building_block_vertices,
        edges,
        lattice_constants,
    ):
        """
        Initialize a :class:`._GraphState` instance.

        Parameters
        ----------
        building_block_vertices : :class:`dict`
            Maps each :class:`.BuildingBlock` to be placed, to a
            :class:`tuple` of :class:`.Vertex` instances, on which
            it should be placed.

        edges : :class:`tuple` of :class:`.Edge`
            The edges which make up the topology graph.

        lattice_constants : :class:`tuple` of :class:`numpy.ndarray`
            A :class:`numpy.ndarray` for each lattice constant.
            Can be an empty :class:`tuple` if the topology graph is
            not periodic.

        """

        self._vertex_building_blocks = {
            vertex.get_id(): building_block
            for building_block, vertices
            in building_block_vertices.items()
            for vertex in vertices
        }
        self._num_building_blocks = {
            building_block: len(vertices)
            for building_block, vertices
            in building_block_vertices.items()
        }
        self._vertices = {
            vertex.get_id(): vertex
            for vertices in building_block_vertices.values()
            for vertex in vertices
        }
        self._edges = edges
        self._lattice_constants = tuple(map(
            np.array,
            lattice_constants,
        ))
        self._vertex_edges = self._get_vertex_edges()

    def _get_vertex_edges(self):
        """
        Get the edges connected to each vertex.

        Returns
        -------
        :class:`dict`
            Maps the id of every vertex to a :class:`tuple` of
            :class:`.Edge` instances connected to it.

        """

        vertex_edges = defaultdict(list)
        for edge in self._edges:
            if edge.is_periodic():
                for vertex_id in edge.get_vertex_ids():
                    periodic_edge = self._get_periodic_edge(
                        edge=edge,
                        reference=vertex_id,
                    )
                    vertex_edges[vertex_id].append(periodic_edge)
            else:
                for vertex_id in edge.get_vertex_ids():
                    vertex_edges[vertex_id].append(edge)
        return vertex_edges

    def _get_periodic_edge(self, edge, reference):
        """
        Get an :class:`.Edge`, with its position correctly set.

        For a periodic edge, its correct position is not at the
        midpoint of the two vertices it connects. Instead, its
        correction position is different for each vertex. The get the
        correct position from the perspective of *vertex1*, *vertex2*
        must first be shifted to its periodic position, and only then
        can the midpoint of the vertices be used to get the edge
        position. An analogous calculation must be done to get the
        position of the edge from the perspective of *vertex2*.

        Parameters
        ----------
        edge : :class:`.Edge`
            The edge whose periodic position must be set.

        reference : :class:`int`
            The id of the vertex, relative to which the edge position
            is being calculated.

        Returns
        -------
        :class:`.Edge`
            A clone of `edge`, shifted to have the correct periodic
            position relative to `reference`.

        """

        vertex1 = self._vertices[reference]
        id1, id2 = edge.get_vertex_ids()
        vertex2 = self._vertices[id1 if reference == id2 else id2]

        direction = 1 if reference == id1 else -1
        periodicity = np.array(edge.get_periodicity())
        end_cell = vertex1.get_cell() + direction*periodicity
        cell_shift = end_cell - vertex2.get_cell()

        shift = sum(
            axis_shift*constant
            for axis_shift, constant in zip(
                cell_shift,
                self._lattice_constants,
            )
        )
        position = (
            (vertex2.get_position()+shift+vertex1.get_position()) / 2
        )
        return edge.with_position(position)

    def clone(self):
        """
        Get a clone.

        Returns
        -------
        :class:`._GraphState`
            The clone. Has the same type as the original instance.

        """

        clone = self.__class__.__new__(self.__class__)
        clone._vertex_building_blocks = dict(
            self._vertex_building_blocks
        )
        clone._vertices = dict(self._vertices)
        clone._vertex_edges = dict(self._vertex_edges)
        clone._edges = self._edges
        clone._lattice_constants = tuple(map(
            np.array,
            self._lattice_constants,
        ))
        clone._num_building_blocks = dict(self._num_building_blocks)
        return clone

    def get_building_block(self, vertex_id):
        """
        Get the building block to be placed on a given vertex.

        Parameters
        ----------
        vertex_id : :class:`int`
            The id of the vertex, on which the building block is to
            be placed.

        Returns
        -------
        :class:`.BuildingBlock`
            The building block.

        """

        return self._vertex_building_blocks[vertex_id]

    def get_vertices(self, vertex_ids=None):
        """
        Yield the topology graph vertices.

        Parameters
        ----------
        vertex_ids : :class:`iterable` of :class:`int`, optional
            The ids of vertices to yield. If ``None``, all vertices
            will be yielded. Can be a single :class:`int` if a
            single vertex is to be yielded.

        Yields
        ------
        :class:`.Vertex`
            A vertex.

        """

        if vertex_ids is None:
            vertex_ids = range(len(self._vertices))
        elif isinstance(vertex_ids, int):
            vertex_ids = (vertex_ids, )

        for vertex_id in vertex_ids:
            yield self._vertices[vertex_id]

    def get_num_vertices(self):
        """
        Get the number of vertices in the topology graph.

        Returns
        -------
        :class:`int`
            The number of vertices in the topology graph.

        """

        return len(self._vertices)

    def get_edge(self, edge_id):
        """
        Get an edge.

        Parameters
        ----------
        edge_id : :class:`int`
            The id of an edge.

        Returns
        -------
        :class:`.Edge`
            An edge.

        """

        return self._edges[edge_id]

    def get_num_edges(self):
        """
        Get the number of edges in the topology graph.

        Returns
        -------
        :class:`int`
            The number of edges.

        """

        return len(self._edges)

    def get_lattice_constants(self):
        """
        Get the lattice constants of the state.

        Returns
        -------
        :class:`tuple` of :class:`numpy.ndarray`
            The lattice constants.

        """

        return tuple(map(np.array, self._lattice_constants))

    def get_edges(self, vertex_id):
        """
        Get the edges connect to a vertex.

        Parameters
        ----------
        vertex_id : class:`int`
            The id of a vertex.

        Returns
        -------
        :class:`tuple` of :class:`.Edge`
            The connected edges.

        """

        return self._vertex_edges[vertex_id]

    def _with_vertices(self, vertices):
        """
        Modify the instance.

        """

        self._vertices = {
            vertex.get_id(): vertex for vertex in vertices
        }
        return self

    def with_vertices(self, vertices):
        """
        Returns a clone holding `vertices`.

        Parameters
        ----------
        vertices : :class:`iterable` of :class:`.Vertex`
            The vertices the clone should hold.

        Returns
        -------
        :class:`._GraphState`
            The clone. Has the same type as the original instance.

        """

        return self.clone()._with_vertices(vertices)

    def _with_lattice_constants(self, lattice_constants):
        """
        Modify the instance.

        """

        self._lattice_constants = tuple(map(
            np.array,
            lattice_constants,
        ))
        return self

    def with_lattice_constants(self, lattice_constants):
        """
        Return a clone holding the `lattice_constants`.

        Parameters
        ----------
        lattice_constants : :class:`tuple` of :class:`numpy.ndarray`
            The lattice constants of the clone. Requires 3 arrays of
            size``(3, )``.

        Returns
        -------
        :class:`._GraphState`
            The clone holding the new lattice constants. Has the same
            type as the original instance.

        """

        return self.clone()._with_lattice_constants(lattice_constants)

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

        return self._num_building_blocks[building_block]

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

        yielded = set()
        for vertex_id in range(max(self._vertex_building_blocks)+1):
            building_block = self._vertex_building_blocks[vertex_id]
            if building_block not in yielded:
                yielded.add(building_block)
                yield building_block
