"""
Covalent Organic Framework
==========================

For usage examples see :class:`.COF`.

#. :class:`Honeycomb`
#. :class:`Hexagonal`
#. :class:`Square`
#. :class:`Kagome`
#. :class:`LinkerlessHoneycomb`

"""

import numpy as np
import itertools as it

from .topology_graph import TopologyGraph, VertexData, Vertex, EdgeData
from ...utilities import vector_angle, flatten


class _COFVertexData(VertexData):
    """
    Holds data for a COF vertex.

    Attributes
    ----------
    id : :class:`int`
        The id of the vertex. Must match the index in
        :attr:`TopologyGraph.vertices`.

    position : :class:`numpy.ndarray`
        The position of the vertex.

    edges : :class:`list` of :class:`.EdgeData`
        The edges connected to the vertex.

    cell : :class:`numpy.ndarray`
        The unit cell in which the vertex is found.

    aligner_edge : :class:`int`
        The edge which is used to align the :class:`.BuildingBlock`
        placed on the vertex. The first :class:`.FunctionalGroup`
        in :attr:`.BuildingBlock.func_groups` is rotated such that
        it lies exactly on this :class:`.Edge`. Must be between
        ``0`` and the number of edges the vertex is connected to.

    """

    def __init__(self, x, y, z):
        """
        Initialize a :class:`.VertexData` instance.

        Parameters
        ----------
        x : :class:`float`
            The x coordinate.

        y : :class:`float`
            The y coordinate.

        z : :class:`float`
            The z coordinate.

        """

        self.aligner_edge = None
        super().__init__(x, y, z)

    @classmethod
    def init_at_center(cls, *vertex_data):
        obj = super().init_at_center(*vertex_data)
        obj.aligner_edge = None
        return obj

    @classmethod
    def init_at_shifted_center(
        cls,
        vertex_data,
        shifts,
        lattice_constants,
        aligner_edge=None
    ):
        """
        Initialize at the center of shifted `vertex_data`.

        Parameters
        ----------
        vertex_data : :class:`tuple` of :class:`.VertexData`
            The vertices at whose center this vertex should be
            intialized.

        shifts : :class:`tuple`
            For every vertex in `vertices` the amount by which it is
            shifted along each axis. For example

            .. code-block:: python

                shifts = (
                    (1, 0, -1),
                    (0, 0, 0)
                )

            means that the first vertex in `vertices` is shifted
            up along the x axis, is not shifted along the y axis
            and is shifted down along the z axis and the second
            vertex is not shifted at all.

        lattice_constants : :class:`tuple` of :class:`numpy.ndarray`
            The a, b and c lattice constants, each written as a vector.

        """

        positions = []
        for vertex, shift in zip(vertex_data, shifts):
            total_shift = 0
            for dim_shift, constant in zip(shift, lattice_constants):
                total_shift += dim_shift * constant
            positions.append(vertex.position + total_shift)

        position = np.divide(
            np.sum(positions, axis=0),
            len(positions)
        )
        return cls(*position)

    def clone(self, clear_edges=False):
        clone = super().clone(clear_edges)
        clone.aligner_edge = self.aligner_edge
        return clone

    def get_vertex(self):
        return _COFVertex(self)


class _COFVertex(Vertex):
    """
    Represents a vertex of a :class:`.COF`.

    Attributes
    ----------
    id : :class:`int`
        The id of the vertex. This should be its index in
        :attr:`.TopologyGraph.vertices`.

    """

    def __init__(self, data):
        """
        Initialize a :class:`_COFVertex`.

        Parameters
        ----------
        data : :class:`_COFVertexData`
            The vertex data.

        """

        # The edge which is used to align the :class:`.BuildingBlock`
        # placed on the vertex. The first :class:`.FunctionalGroup`
        # in :attr:`.BuildingBlock.func_groups` is rotated such that
        # it lies exactly on this :class:`.Edge`. Must be between
        # ``0`` and the number of edges the vertex is connected to.
        self._aligner_edge = data.aligner_edge
        super().__init__(data)

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
        clone._aligner_edge = self._aligner_edge
        return clone

    def get_aligner_edge(self):
        return self._aligner_edge

    def place_building_block(self, building_block, vertices, edges):
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

        if len(building_block.func_groups) == 2:
            return self._place_linear_building_block(
                building_block=building_block,
                vertices=vertices,
                edges=edges
            )
        return self._place_nonlinear_building_block(
            building_block=building_block,
            vertices=vertices,
            edges=edges
        )

    def _place_linear_building_block(
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
        fg_centroid = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_bonder_ids()
        )
        start = fg_centroid - self._position
        e0_coord = (
            edges[self._edge_ids[0]].get_position(self, vertices)
        )
        e1_coord = (
            edges[self._edge_ids[1]].get_position(self, vertices)
        )
        target = e0_coord - e1_coord

        if self._edge_ids[self._aligner_edge] != self._edge_ids[0]:
            target *= -1

        building_block.apply_rotation_between_vectors(
            start=start,
            target=target,
            origin=self._position
        )
        start = building_block.get_centroid_centroid_direction_vector()
        building_block.apply_rotation_to_minimize_angle(
            start=start,
            target=self._position,
            axis=target,
            origin=self._position,
        )
        return building_block.get_position_matrix()

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
            target=[0, 0, 1],
            origin=self._position
        )
        fg_bonder_centroid = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_bonder_ids()
        )
        start = fg_bonder_centroid - self._position

        aligner_edge = edges[self._edge_ids[self._aligner_edge]]
        edge_coord = aligner_edge.get_position(self, vertices)
        target = edge_coord - self._position
        building_block.apply_rotation_to_minimize_angle(
            start=start,
            target=target,
            axis=[0, 0, 1],
            origin=self._position
        )
        return building_block.get_position_matrix()

    def assign_func_groups_to_edges(
        self,
        building_block,
        vertices,
        edges
    ):
        """
        Assign functional groups to edges.

        Each :class:`.FunctionalGroup` of the `building_block` needs
        to be associated with one of the :class:`.Edge` instances in
        :attr:`edges`.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is needs to have
            functional groups assigned to edges.

        vertices : :class:`tuple` of :class:`.Vertex`
            All vertices in the topology graph. The index of each
            vertex must match its :class:`~.Vertex.id`.

        edges : :class:`tuple` of :class:`.Edge`
            All edges in the topology graph. The index of each
            edge must match its :class:`~.Edge.id`.

        Returns
        -------
        :class:`dict`
            A mapping from the id of a functional group in
            `building_block` to the id of the edge in :attr:`edges` it
            is assigned to.

        """

        if len(building_block.func_groups) == 2:
            return self._assign_func_groups_to_linear_edges(
                building_block=building_block,
                vertices=vertices,
                edges=edges
            )
        return self._assign_func_groups_to_nonlinear_edges(
                building_block=building_block,
                vertices=vertices,
                edges=edges
            )

    def _assign_func_groups_to_linear_edges(
        self,
        building_block,
        vertices,
        edges
    ):
        return {
            fg_id: edge_id for fg_id, edge_id in enumerate(sorted(
                self._edge_ids,
                key=self._get_fg0_distance(
                    building_block=building_block,
                    vertices=vertices,
                    edges=edges
                )
            ))
        }

    def _get_fg0_distance(self, building_block, vertices, edges):
        fg_coord = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_bonder_ids()
        )

        def distance(edge_id):
            displacement = edges[edge_id].get_position(
                reference=self,
                vertices=vertices
            ) - fg_coord
            return np.linalg.norm(displacement)

        return distance

    def _assign_func_groups_to_nonlinear_edges(
        self,
        building_block,
        vertices,
        edges
    ):
        # The idea is to order the functional groups in building_block
        # by their angle from func_groups[0] and the bonder centroid,
        #  going in the clockwise direction.
        #
        # The edges are also ordered by their angle from aligner_edge
        # and the edge centroid going in the clockwise direction.
        #
        # Once the fgs and edges are ordered, zip and assign them.

        fg0_coord = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_bonder_ids()
        )
        bonder_centroid = building_block.get_centroid(
            atom_ids=building_block.get_bonder_ids()
        )
        fg0_direction = fg0_coord-bonder_centroid
        axis = np.cross(
            fg0_direction,
            building_block.get_bonder_plane_normal()
        )
        func_groups = sorted(
            range(len(building_block.func_groups)),
            key=self._get_func_group_angle(
                building_block=building_block,
                fg0_direction=fg0_direction,
                bonder_centroid=bonder_centroid,
                axis=axis
            )
        )
        edge_ids = sorted(
            self._edge_ids,
            key=self._get_edge_angle(axis, vertices, edges)
        )
        assignments = {}
        for edge_id, fg_id in zip(edge_ids, func_groups):
            assignments[fg_id] = edge_id
        return assignments

    @staticmethod
    def _get_func_group_angle(
        building_block,
        fg0_direction,
        bonder_centroid,
        axis
    ):

        def angle(fg_id):
            func_group = building_block.func_groups[fg_id]
            coord = building_block.get_centroid(
                atom_ids=func_group.get_bonder_ids()
            )
            fg_direction = coord-bonder_centroid
            theta = vector_angle(fg0_direction, fg_direction)

            projection = fg_direction @ axis
            if theta > 0 and projection < 0:
                return 2*np.pi - theta
            return theta

        return angle

    def _get_edge_angle(self, axis, vertices, edges):

        aligner_edge = edges[self._edge_ids[self._aligner_edge]]
        aligner_edge_coord = aligner_edge.get_position(self, vertices)
        connected_edges = tuple(edges[id_] for id_ in self._edge_ids)
        edge_centroid = self._get_edge_centroid(
            centroid_edges=connected_edges,
            vertices=vertices
        )
        # This axis is used to figure out the clockwise direction.
        aligner_edge_direction = aligner_edge_coord - edge_centroid

        def angle(edge_id):
            coord = edges[edge_id].get_position(self, vertices)
            edge_direction = coord - edge_centroid
            theta = vector_angle(
                vector1=edge_direction,
                vector2=aligner_edge_direction
            )

            projection = edge_direction @ axis
            if theta > 0 and projection < 0:
                return 2*np.pi - theta
            return theta

        return angle

    def __str__(self):
        return (
            f'Vertex(id={self.id}, '
            f'position={self._position.tolist()}, '
            f'aligner_edge={self._aligner_edge})'
        )


class COF(TopologyGraph):
    """
    Represents a COF topology graph.

    COF topologies are added by creating a subclass which defines the
    :attr:`vertices` and :attr:`edges` of the topology as class
    attributes.

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

    Examples
    --------
    :class:`COF` instances can be made by supplying only
    the lattice size (using :class:`.Honeycomb` as an example)

    .. code-block:: python

        import stk

        bb1 = stk.BuildingBlock('NCCN', ['amine'])
        bb2 = stk.BuildingBlock('O=CC(C=O)C=O', ['aldehyde'])
        cof1 = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topology_graph=stk.cof.Honeycomb((2, 2, 1))
        )

    Different structural isomers of COFs can be made by using the
    `vertex_alignments` optional parameter

    .. code-block:: python

        lattice = stk.cof.Honeycomb(
            lattice_size=(2, 2, 1),
            vertex_alignments={0: 1, 2: 2}
        )
        cof2 = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topology_graph=lattice
        )

    The parameter maps the :attr:`~.Vertex.id` of a vertex to a number
    between 0 (inclusive) and the number of edges the vertex is
    connected to (exclusive). So a vertex connected to three edges
    can be mapped to ``0``, ``1`` or ``2``.

    You can also build COFs with multiple building blocks, but you
    have to assign each building block to a vertex with
    `building_block_vertices`.

    .. code-block:: python

        lattice = stk.cof.Honeycomb(
            lattice_size=(2, 2, 1),
            vertex_alignments={0: 1, 2: 2}
        )
        bb3 = stk.BuildingBlock('NCOCN', ['amine'])
        cof2 = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2, bb3],
            topology_graph=lattice
            building_block_vertices={
                bb1: lattice.verices[:2],
                bb2: lattice.verices[4:],
                bb3: lattice.verices[2:4]
            }
        )

    """

    def __init_subclass__(cls, **kwargs):
        for i, vertex in enumerate(cls.vertex_data):
            vertex.id = i
        for i, edge in enumerate(cls.edge_data):
            edge.id = i
            edge.lattice_constants = tuple(
                np.array(constant)
                for constant in cls._lattice_constants
            )
        return super().__init_subclass__(**kwargs)

    def __init__(
        self,
        lattice_size,
        periodic=False,
        vertex_alignments=None,
        num_processes=1
    ):
        """
        Initialize a :class:`.COF`.

        Parameters
        ----------
        lattice_size : :class:`tuple` of :class:`int`
            The number of unit cells which should be placed along the
            x, y and z dimensions, respectively.

        periodic : :class:`bool`, optional
            If periodic bonds are to be made across the lattice,
            this should be ``True``. If ``False`` the functional
            groups on the ends of the lattice will be unreacted.

        vertex_alignments : :class:`dict`, optional
            A mapping from the :attr:`.Vertex.id` of a :class:`.Vertex`
            :attr:`vertices` to an :class:`.Edge` connected to it.
            The :class:`.Edge` is used to align the first
            :class:`.FunctionalGroup` of a :class:`.BuildingBlock`
            placed on that vertex. Only vertices which need to have
            their default edge changed need to be present in the
            :class:`dict`. If ``None`` then the default edge is used
            for each vertex. Changing which :class:`.Edge` is used will
            mean that the topology graph represents different
            structural isomers. The edge is refered to by a number
            between ``0`` (inclusive) and the number of edges the
            vertex is connected to (exclusive).

        num_processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        """

        if vertex_alignments is None:
            vertex_alignments = {}

        self._lattice_size = lattice_size
        self._periodic = periodic

        vertex_data = self._get_vertex_data(vertex_alignments)
        edge_data = self._get_edge_data(vertex_data)

        vertex_data = tuple(
            vertex
            for clones in flatten(vertex_data, {dict})
            for vertex in clones.values()
        )
        super().__init__(vertex_data, edge_data, (), num_processes)

    def _get_vertex_data(self, vertex_alignments):
        """
        Create the vertex data of the topology graph instance.

        Parameters
        ---------
        vertex_alignments : :class:`dict`
            A mapping from the :attr:`.Vertex.id` of a :class:`.Vertex`
            :attr:`vertices` to an :class:`.Edge` connected to it.
            The :class:`.Edge` is used to align the first
            :class:`.FunctionalGroup` of a :class:`.BuildingBlock`
            placed on that vertex. Only vertices which need to have
            their default edge changed need to be present in the
            :class:`dict`. If ``None`` then the default edge is used
            for each vertex. Changing which :class:`.Edge` is used will
            mean that the topology graph represents different
            structural isomers. The edge is refered to by a number
            between ``0`` (inclusive) and the number of edges the
            vertex is connected to (exclusive).

        Returns
        -------
        :class:`list`
            A nested :class:`list` which can be indexed as
            ``vertices[x][y][z]``, which will return a :class:`dict`
            for the unit cell at (x, y, z). The :class:`dict` maps
            the vertices in :attr:`vertex_data` to its clone for that
            unit cell.

        """

        xdim, ydim, zdim = (range(dim) for dim in self._lattice_size)
        # vertex_clones is indexed as vertex_clones[x][y][z]
        vertex_clones = [
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
        vertices = it.product(cells, self.vertex_data)
        for cell, vertex in vertices:
            x, y, z = cell
            clone = vertex.clone(True)
            clone.cell = np.array(cell)
            clone.aligner_edge = vertex_alignments.get(vertex.id, 0)
            # Shift the clone so that it's within the cell.
            for axis, dim in zip(cell, self._lattice_constants):
                clone.position += axis * dim

            vertex_clones[x][y][z][vertex] = clone
        return vertex_clones

    def _get_edge_data(self, vertex_data):
        """
        Create the edge data of the topology graph instance.

        Parameters
        ----------
        vertex_data : :class:`list`
            A nested :class:`list` which can be indexed as
            ``vertex_data[x][y][z]``, which will return a :class:`dict`
            for the unit cell at (x, y, z). The :class:`dict` maps
            the vertices in :attr:`vertex_data` to the clones for that
            unit cell.

        Returns
        -------
        :class:`tuple` of :class:`.EdgeData`
            The edge data of the topology graph instance.

        """

        edge_clones = []
        # Make a clone for each edge for each unit cell.
        xdim, ydim, zdim = (range(dim) for dim in self._lattice_size)
        cells = it.product(xdim, ydim, zdim)
        edges = it.product(cells, self.edge_data)
        for cell, edge in edges:
            x, y, z = cell
            # The cell in which the second vertex of the edge is found.
            periodic_cell = np.array(cell) + edge.periodicity
            # Wrap around periodic cells, ie those that are less than 0
            # or greater than the lattice size along any dimension.
            dims = zip(periodic_cell, self._lattice_size)
            cell2_x, cell2_y, cell2_z = np.array([
                (dim+max_dim) % max_dim
                for dim, max_dim in dims
            ])
            # Make a vertex map which accounts for the fact that
            # v1 is in cell2.
            v0, v1 = edge.vertices
            vertex_map = {
                v0: vertex_data[x][y][z][v0],
                v1: vertex_data[cell2_x][cell2_y][cell2_z][v1]
            }
            # If the edge is not periodic if periodic_cell is did not
            # have to wrap around.
            dims = zip(periodic_cell, self._lattice_size)
            edge_is_not_periodic = all(
                dim >= 0 and dim < max_dim
                for dim, max_dim in dims
            )
            clone = edge.clone(vertex_map, True, True)
            edge_clones.append(clone)
            if edge_is_not_periodic:
                clone.periodicity = np.array([0, 0, 0])

        return tuple(edge_clones)

    def _before_react(self, mol, vertices, edges):
        if self._periodic:
            return vertices, edges
        return vertices, tuple(
            edge for edge in edges if not edge.is_periodic()
        )

    def assign_building_blocks_to_vertices(self, building_blocks):
        """
        Assign `building_blocks` to :attr:`vertices`.

        Parameters
        ----------
        building_blocks : :class:`list` of :class:`.Molecule`
            The :class:`.BuildingBlock` and
            :class:`ConstructedMolecule` instances which
            represent the building block molecules used for
            construction. Only one instance is present per building
            block molecule, even if multiples of that building block
            join up to form the :class:`ConstructedMolecule`.

        Returns
        -------
        :class:`dict`
            Maps the `building_blocks`, to the
            :class:`~.topologies.base.Vertex` objects in
            :attr:`vertices` they are placed on during construction.
            The :class:`dict` has the form

            .. code-block:: python

                building_block_vertices = {
                    BuildingBlock(...): [Vertex(...), Vertex(...)],
                    BuildingBlock(...): [
                        Vertex(...),
                        Vertex(...),
                        Vertex(...),
                    ]
                    ConstructedMolecule(...): [Vertex(...)]
                }

        Raises
        ------
        :class:`ValueError`
            If there is more than one building with a given number
            of functional groups.

        """

        bb_by_degree = {}
        for bb in building_blocks:
            num_fgs = len(bb.func_groups)
            if num_fgs in bb_by_degree:
                raise ValueError(
                    'If there are multiple building blocks with the '
                    'same number of functional groups, '
                    'building_block_vertices must be set explicitly.'
                )
            bb_by_degree[num_fgs] = bb

        building_block_vertices = {}
        for vertex in self.vertices:
            bb = bb_by_degree[vertex.get_num_edges()]
            building_block_vertices[bb] = (
                building_block_vertices.get(bb, [])
            )
            building_block_vertices[bb].append(vertex)
        return building_block_vertices

    def _get_scale(self, mol):
        """
        Get the scale used for the positions of :attr:`vertices`.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The molecule being constructed.

        Returns
        -------
        :class:`float` or :class:`list` of :class:`float`
            The value by which the position of each :class:`Vertex` is
            scaled. Can be a single number if all axes are scaled by
            the same amount or a :class:`list` of three numbers if
            each axis is scaled by a different value.

        """

        return 5*max(
            bb.get_maximum_diameter()
            for bb in mol.building_block_vertices
        )

    def __repr__(self):
        vertex_alignments = ', '.join(
            f'{v.id}: {v.get_aligner_edge()}'
            # Only get the vertices in the first unit cell.
            for v in self.vertices[:len(self.vertex_data)]
        )

        x, y, z = self._lattice_size
        periodic = ', periodic=True' if self._periodic else ''
        return (
            f'cof.{self.__class__.__name__}('
            f'lattice_size=({x}, {y}, {z}), '
            f'vertex_alignments={{{vertex_alignments}}}'
            f'{periodic})'
        )


class Honeycomb(COF):
    """
    Represents a honeycomb COF topology graph.

    Building blocks with three and two functional groups are required
    for this topology graph.

    See :class:`.COF` for more details and examples.

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

    _lattice_constants = _a, _b, _c = (
        np.array([1., 0., 0.]),
        np.array([0.5, 0.866, 0]),
        np.array([0, 0, 5/1.7321])
    )

    _vertex_data = (
        _COFVertexData(*((1/3)*_a + (1/3)*_b + (1/2)*_c)),
        _COFVertexData(*((2/3)*_a + (2/3)*_b + (1/2)*_c))
    )

    vertex_data = (
        *_vertex_data,
        _COFVertexData.init_at_center(
            _vertex_data[0], _vertex_data[1]
        ),
        _COFVertexData.init_at_shifted_center(
            vertex_data=(_vertex_data[0], _vertex_data[1]),
            shifts=((0, 0, 0), (0, -1, 0)),
            lattice_constants=_lattice_constants
        ),
        _COFVertexData.init_at_shifted_center(
            vertex_data=(_vertex_data[0], _vertex_data[1]),
            shifts=((0, 0, 0), (-1, 0, 0)),
            lattice_constants=_lattice_constants
        )
    )

    edge_data = (
        EdgeData(vertex_data[2], vertex_data[0]),
        EdgeData(vertex_data[2], vertex_data[1]),

        EdgeData(vertex_data[3], vertex_data[0]),
        EdgeData(
            vertex_data[3],
            vertex_data[1],
            periodicity=(0, -1, 0)
        ),

        EdgeData(vertex_data[4], vertex_data[0]),
        EdgeData(
            vertex_data[4],
            vertex_data[1],
            periodicity=(-1, 0, 0)
        )
    )


class Hexagonal(COF):
    """
    Represents a hexagonal COF topology graph.

    Building blocks with six and two functional groups are required
    for this topology graph.

    See :class:`.COF` for more details and examples.

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

    _lattice_constants = _a, _b, _c = (
        np.array([1., 0., 0.]),
        np.array([0.5, 0.866, 0]),
        np.array([0, 0, 5/1.7321])
    )

    _vertex_data = (
        _COFVertexData(*((1/4)*_a + (1/4)*_b + (1/2)*_c)),
        _COFVertexData(*((1/4)*_a + (3/4)*_b + (1/2)*_c)),
        _COFVertexData(*((3/4)*_a + (1/4)*_b + (1/2)*_c)),
        _COFVertexData(*((3/4)*_a + (3/4)*_b + (1/2)*_c))
    )

    vertex_data = (
        *_vertex_data,
        _COFVertexData.init_at_center(
            _vertex_data[0], _vertex_data[1]
        ),
        _COFVertexData.init_at_center(
            _vertex_data[0], _vertex_data[2]
        ),
        _COFVertexData.init_at_center(
            _vertex_data[1], _vertex_data[2]
        ),
        _COFVertexData.init_at_center(
            _vertex_data[1], _vertex_data[3]
        ),
        _COFVertexData.init_at_center(
            _vertex_data[2], _vertex_data[3]
        ),
        _COFVertexData.init_at_shifted_center(
            vertex_data=(_vertex_data[0], _vertex_data[2]),
            shifts=((0, 0, 0), (-1, 0, 0)),
            lattice_constants=_lattice_constants
        ),
        _COFVertexData.init_at_shifted_center(
            vertex_data=(_vertex_data[0], _vertex_data[1]),
            shifts=((0, 0, 0), (0, -1, 0)),
            lattice_constants=_lattice_constants
        ),
        _COFVertexData.init_at_shifted_center(
            vertex_data=(_vertex_data[0], _vertex_data[3]),
            shifts=((0, 0, 0), (0, -1, 0)),
            lattice_constants=_lattice_constants
        ),
        _COFVertexData.init_at_shifted_center(
            vertex_data=(_vertex_data[2], _vertex_data[1]),
            shifts=((0, 0, 0), (1, -1, 0)),
            lattice_constants=_lattice_constants
        ),
        _COFVertexData.init_at_shifted_center(
            vertex_data=(_vertex_data[2], _vertex_data[3]),
            shifts=((0, 0, 0), (0, -1, 0)),
            lattice_constants=_lattice_constants
        ),
        _COFVertexData.init_at_shifted_center(
            vertex_data=(_vertex_data[1], _vertex_data[3]),
            shifts=((0, 0, 0), (-1, 0, 0)),
            lattice_constants=_lattice_constants
        ),
        _COFVertexData.init_at_shifted_center(
            vertex_data=(_vertex_data[3], _vertex_data[0]),
            shifts=((0, 0, 0), (1, 0, 0)),
            lattice_constants=_lattice_constants
        )
    )

    edge_data = (
        EdgeData(vertex_data[4], vertex_data[0]),
        EdgeData(vertex_data[4], vertex_data[1]),

        EdgeData(vertex_data[5], vertex_data[0]),
        EdgeData(vertex_data[5], vertex_data[2]),

        EdgeData(vertex_data[6], vertex_data[1]),
        EdgeData(vertex_data[6], vertex_data[2]),

        EdgeData(vertex_data[7], vertex_data[1]),
        EdgeData(vertex_data[7], vertex_data[3]),

        EdgeData(vertex_data[8], vertex_data[2]),
        EdgeData(vertex_data[8], vertex_data[3]),

        EdgeData(vertex_data[9], vertex_data[0]),
        EdgeData(
            vertex_data[9],
            vertex_data[2],
            periodicity=(-1, 0, 0)
        ),

        EdgeData(vertex_data[10], vertex_data[0]),
        EdgeData(
            vertex_data[10],
            vertex_data[1],
            periodicity=(0, -1, 0)
        ),

        EdgeData(vertex_data[11], vertex_data[0]),
        EdgeData(
            vertex_data[11],
            vertex_data[3],
            periodicity=(0, -1, 0)
        ),

        EdgeData(vertex_data[12], vertex_data[2]),
        EdgeData(
            vertex_data[12],
            vertex_data[1],
            periodicity=(1, -1, 0)
        ),

        EdgeData(vertex_data[13], vertex_data[2]),
        EdgeData(
            vertex_data[13],
            vertex_data[3],
            periodicity=(0, -1, 0)
        ),

        EdgeData(vertex_data[14], vertex_data[1]),
        EdgeData(
            vertex_data[14],
            vertex_data[3],
            periodicity=(-1, 0, 0)
        ),

        EdgeData(vertex_data[15], vertex_data[3]),
        EdgeData(
            vertex_data[15],
            vertex_data[0],
            periodicity=(1, 0, 0)
        )
    )


class Square(COF):
    """
    Represents a sqaure COF topology graph.

    Building blocks with four and two functional groups are required
    for this topology graph.

    See :class:`.COF` for more details and examples.

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

    _lattice_constants = _a, _b, _c = (
        np.array([1., 0., 0.]),
        np.array([0., 1., 0.]),
        np.array([0., 0., 1.])
    )

    _vertex_data = (
        _COFVertexData(*((0.5)*_a + (0.5)*_b + (0.5)*_c)),
    )
    vertex_data = (
        *_vertex_data,
        _COFVertexData.init_at_shifted_center(
            vertex_data=(_vertex_data[0], _vertex_data[0]),
            shifts=((0, 0, 0), (1, 0, 0)),
            lattice_constants=_lattice_constants
        ),
        _COFVertexData.init_at_shifted_center(
            vertex_data=(_vertex_data[0], _vertex_data[0]),
            shifts=((0, 0, 0), (0, 1, 0)),
            lattice_constants=_lattice_constants
        )

    )

    edge_data = (
        EdgeData(vertex_data[1], vertex_data[0]),
        EdgeData(
            vertex_data[1],
            vertex_data[0],
            periodicity=(1, 0, 0)
        ),
        EdgeData(vertex_data[2], vertex_data[0]),
        EdgeData(
            vertex_data[2],
            vertex_data[0],
            periodicity=(0, 1, 0)
        )
    )


class Kagome(COF):
    """
    Represents a kagome COF topology graph.

    Building blocks with four and two functional groups are required
    for this topology graph.

    See :class:`.COF` for more details and examples.

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

    _lattice_constants = _a, _b, _c = (
        np.array([1., 0., 0.]),
        np.array([0.5, 0.866, 0.]),
        np.array([0., 0., 5/1.7321])
    )

    _vertex_data = (
        _COFVertexData(*((1/4)*_a + (3/4)*_b + (0.5)*_c)),
        _COFVertexData(*((3/4)*_a + (3/4)*_b + (1/2)*_c)),
        _COFVertexData(*((3/4)*_a + (1/4)*_b + (1/2)*_c))
    )

    vertex_data = (
        *_vertex_data,
        _COFVertexData.init_at_center(
            _vertex_data[0],
            _vertex_data[1]
        ),
        _COFVertexData.init_at_center(
            _vertex_data[0],
            _vertex_data[2]
        ),
        _COFVertexData.init_at_center(
            _vertex_data[1],
            _vertex_data[2]
        ),
        _COFVertexData.init_at_shifted_center(
            vertex_data=(_vertex_data[0], _vertex_data[1]),
            shifts=((0, 0, 0), (-1, 0, 0)),
            lattice_constants=_lattice_constants
        ),
        _COFVertexData.init_at_shifted_center(
            vertex_data=(_vertex_data[0], _vertex_data[2]),
            shifts=((0, 0, 0), (-1, 1, 0)),
            lattice_constants=_lattice_constants
        ),
        _COFVertexData.init_at_shifted_center(
            vertex_data=(_vertex_data[1], _vertex_data[2]),
            shifts=((0, 0, 0), (0, 1, 0)),
            lattice_constants=_lattice_constants
        )

    )

    edge_data = (
        EdgeData(vertex_data[3], vertex_data[0]),
        EdgeData(vertex_data[3], vertex_data[1]),

        EdgeData(vertex_data[4], vertex_data[0]),
        EdgeData(vertex_data[4], vertex_data[2]),

        EdgeData(vertex_data[5], vertex_data[1]),
        EdgeData(vertex_data[5], vertex_data[2]),

        EdgeData(vertex_data[6], vertex_data[0]),
        EdgeData(
            vertex_data[6],
            vertex_data[1],
            periodicity=(-1, 0, 0)
        ),

        EdgeData(vertex_data[7], vertex_data[0]),
        EdgeData(
            vertex_data[7],
            vertex_data[2],
            periodicity=(-1, 1, 0)
        ),

        EdgeData(vertex_data[8], vertex_data[1]),
        EdgeData(
            vertex_data[8],
            vertex_data[2],
            periodicity=(0, 1, 0)
        )
    )


class LinkerlessHoneycomb(COF):
    """
    Represents a honeycomb COF topology graph.

    Building blocks with three functional groups are required
    for this topology graph.

    See :class:`.COF` for more details and examples.

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

    _lattice_constants = _a, _b, _c = (
        np.array([1., 0., 0.]),
        np.array([0.5, 0.866, 0.]),
        np.array([0., 0., 5/1.7321])
    )

    vertex_data = (
        _COFVertexData(*((1/3)*_a + (1/3)*_b + (1/2)*_c)),
        _COFVertexData(*((2/3)*_a + (2/3)*_b + (1/2)*_c))
    )

    edge_data = (
        EdgeData(vertex_data[0], vertex_data[1]),
        EdgeData(
            vertex_data[0],
            vertex_data[1],
            periodicity=(-1, 0, 0)
        ),
        EdgeData(
            vertex_data[0],
            vertex_data[1],
            periodicity=(0, -1, 0)
        )
    )
