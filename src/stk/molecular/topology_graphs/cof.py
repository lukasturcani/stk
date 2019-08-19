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

from .topology_graph import TopologyGraph, Vertex, Edge
from ...utilities import vector_angle, flatten


class _COFVertex(Vertex):
    """
    Represents a vertex of a :class:`.COF`.

    Attributes
    ----------
    id : :class:`int`
        The id of the vertex. This should be its index in
        :attr:`TopologyGraph.vertices`.

    edges : :class:`list` of :class:`.Edge`
        The edges the :class:`Vertex` is connected to.

    aligner_edge : :class:`.Edge`
        The :class:`.Edge` in :attr:`edges`, which is used to align the
        :class:`.BuildingBlock` placed on the vertex. The first
        :class:`.FunctionalGroup` in :attr:`.BuildingBlock.func_groups`
        is rotated such that it lies exactly on this :class:`.Edge`.

    """

    def __init__(self, x, y, z):
        """
        Initialize a :class:`_COFVertex`.

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
    def init_at_center(cls, *vertices):
        """
        Initialize at the center of `vertices`.

        Parameters
        ----------
        vertices : :class:`.Vertex`
            Vertices at whose center this vertex should be initialized.

        Returns
        -------
        :class:`.Vertex`
            The vertex.

        """

        center = sum(vertex.get_position() for vertex in vertices)
        center /= len(vertices)
        return cls(*center)

    @classmethod
    def init_at_shifted_center(
        cls,
        vertices,
        shifts,
        lattice_constants
    ):
        """
        Initialize at the center of shifted `vertices`.

        Parameters
        ----------
        vertices : :class:`tuple` of :class:`_COFVertex`
            The vertics at whose center this vertex should be
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
        for vertex, shift in zip(vertices, shifts):
            total_shift = 0
            for dim_shift, constant in zip(shift, lattice_constants):
                total_shift += dim_shift * constant
            position = vertex.get_position() + total_shift
            positions.append(position)

        position = np.divide(
            np.sum(positions, axis=0),
            len(positions)
        )
        return cls(*position)

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
        if self.aligner_edge is not None:
            clone.aligner_edge = self.aligner_edge.clone(
                add_to_vertices=False
            )
        else:
            clone.aligner_edge = None
        return clone

    def apply_scale(self, scale):
        """
        Scale the position by `scale`.

        Parameters
        ----------
        scale : :class:`float` or :class:`list`of :class:`float`
            The value by which the position of the :class:`Vertex` is
            scaled. Can be a single number if all axes are scaled by
            the same amount or a :class:`list` of three numbers if
            each axis is scaled by a different value.

        Returns
        -------
        :class:`Vertex`
            The vertex is returned.

        """

        self._position *= scale
        self.aligner_edge.apply_scale(scale)
        return self

    def place_building_block(self, building_block):
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

        if len(building_block.func_groups) == 2:
            return self._place_linear_building_block(building_block)
        return self._place_nonlinear_building_block(building_block)

    def _place_linear_building_block(self, building_block):
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
        fg_centroid = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_bonder_ids()
        )
        start = fg_centroid - self._position
        e0_coord = self.edges[0].get_position(self)
        e1_coord = self.edges[1].get_position(self)
        target = e0_coord - e1_coord

        if self.aligner_edge.id != self.edges[0].id:
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
            target=[0, 0, 1],
            origin=self._position
        )
        fg_bonder_centroid = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_bonder_ids()
        )
        start = fg_bonder_centroid - self._position
        edge_coord = self.aligner_edge.get_position(self)
        target = edge_coord - self._position
        building_block.apply_rotation_to_minimize_angle(
            start=start,
            target=target,
            axis=[0, 0, 1],
            origin=self._position
        )
        return building_block.get_position_matrix()

    def assign_func_groups_to_edges(self, building_block):
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

        Returns
        -------
        :class:`dict`
            A mapping from the id of a functional group in
            `building_block` to the id of the edge in :attr:`edges` it
            is assigned to.

        """

        if len(building_block.func_groups) == 2:
            return self._assign_func_groups_to_linear_edges(
                building_block=building_block
            )
        return self._assign_func_groups_to_nonlinear_edges(
                building_block=building_block
            )

    def _assign_func_groups_to_linear_edges(self, building_block):
        return {
            fg_id: e.id for fg_id, e in enumerate(sorted(
                self.edges,
                key=self._get_fg0_distance(building_block)
            ))
        }

    def _get_fg0_distance(self, building_block):
        fg_coord = building_block.get_centroid(
            atom_ids=building_block.func_groups[0].get_bonder_ids()
        )

        def distance(edge):
            displacement = edge.get_position(self) - fg_coord
            return np.linalg.norm(displacement)

        return distance

    def _assign_func_groups_to_nonlinear_edges(self, building_block):
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
        edges = sorted(self.edges, key=self._get_edge_angle(axis))
        assignments = {}
        for edge, fg_id in zip(edges, func_groups):
            assignments[fg_id] = edge.id
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

    def _get_edge_angle(self, axis):

        aligner_edge_coord = self.aligner_edge.get_position(self)
        edge_centroid = self._get_edge_centroid()
        # This axis is used to figure out the clockwise direction.
        aligner_edge_direction = aligner_edge_coord - edge_centroid

        def angle(edge):
            coord = edge.get_position(self)
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
        x, y, z = self._position
        if self.aligner_edge is not None:
            aligner_edge = self.edges.index(self.aligner_edge)
        else:
            aligner_edge = None
        return (
            f'Vertex(id={self.id}, '
            f'position={[x, y, z]}, '
            f'aligner_edge={aligner_edge})'
        )


class COF(TopologyGraph):
    """
    Represents a COF topology graph.

    COF topologies are added by creating a subclass which defines the
    :attr:`vertices` and :attr:`edges` of the topology as class
    attributes.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    Examples
    --------
    :class:`COF` instances can be made without supplying
    additional arguments (using :class:`.Honeycomb` as an example)

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

        v0 = stk.cof.Honeycomb.vertices[0]
        v2 = stk.cof.Honeycomb.vertices[2]
        lattice = stk.cof.Honeycomb(
            lattice_size=(2, 2, 1),
            vertex_alignments={
                v0: v0.edges[1],
                v2: v2.edges[2]
            }
        )
        cof2 = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2],
            topology_graph=lattice
        )

    By changing which edge each vertex is aligned with, a different
    structural isomer of the COF can be formed.

    Note the in the `vertex_alignments` parameter the class vertices
    and edges are used, however when the `building_block_vertices`
    parameter is used, the instance vertices are used. **These are not
    interchangeable!**

    .. code-block:: python

        # Use the class vertices and edges to set vertex_alignments
        # and create a topology graph.
        v0 = stk.cof.Honeycomb.vertices[0]
        v2 = stk.cof.Honeycomb.vertices[2]
        lattice = stk.cof.Honeycomb(
            lattice_size=(2, 2, 1),
            vertex_alignments={
                v0: v0.edges[1],
                v2: v2.edges[2]
            }
        )
        bb3 = stk.BuildingBlock('NCOCN', ['amine'])
        cof2 = stk.ConstructedMolecule(
            building_blocks=[bb1, bb2, bb3],
            topology_graph=lattice
            # Use the instance vertices in the building_block_vertices
            # parameter.
            building_block_vertices={
                bb1: lattice.vertices[:2],
                bb2: lattice.vertices[4:],
                bb3: lattice.vertices[2:4]
            }
        )

    The example above also demonstrates how COFs with many building
    blocks can be built. You can add as many :class:`.BuildingBlock`
    instances into `building_blocks` as you like. If you do not
    assign where each building block is placed with
    `building_block_vertices`, they will be placed on the
    :attr:`vertices` of the :class:`.COF` at random. Random
    placement will account for the fact that the length of
    :attr:`.BuildingBlock.func_groups` needs to match the number of
    edges connected to a vertex.

    """

    def __init_subclass__(cls, **kwargs):
        for i, vertex in enumerate(cls.vertices):
            vertex.id = i
        for i, edge in enumerate(cls.edges):
            edge.id = i
            edge._lattice_constants = tuple(
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
            A mapping from a :class:`.Vertex`, in the class
            attribute :attr:`vertices`, to an :class:`.Edge`, in the
            class attribute :attr:`edges`, connected to it.
            The :class:`.Edge` is used to align the first
            :class:`.FunctionalGroup` of a :class:`.BuildingBlock`
            placed on that vertex. Only vertices which need to have
            their default alignment changed need to be present in the
            :class:`dict`. If ``None`` then the first :class:`.Edge`
            in :attr:`.Vertex.edges` is used for each vertex. Changing
            which :class:`.Edge` is used will mean that the topology
            graph will represent a different structural isomer.

            The vertices and edges can also be referred to by their
            indices. The vertices are referred to by their index in
            the class attribute :attr:`vertices` while the edges are
            referred to by their index in :attr:`.Vertex.edges`.

        num_processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        """

        if vertex_alignments is None:
            vertex_alignments = {}
        vertex_alignments = self._normalize_vertex_alignments(
            vertex_alignments=vertex_alignments
        )

        self._lattice_size = lattice_size
        self._periodic = periodic

        vertices = self._get_instance_vertices(vertex_alignments)
        edges = self._get_instance_edges(vertices)

        vertices = tuple(
            vertex
            for clones in flatten(vertices, {dict})
            for vertex in clones.values()
        )
        super().__init__(vertices, edges, (), num_processes)

    def _normalize_vertex_alignments(self, vertex_alignments):
        """
        Normalize different `vertex_alignments` input forms.

        Parameters
        ----------
        vertex_alignments : :class:`dict`
            A mapping from a :class:`.Vertex`, in the class
            attribute :attr:`vertices`, to an :class:`.Edge`, in the
            class attribute :attr:`edges`, connected to it.
            The :class:`.Edge` is used to align the first
            :class:`.FunctionalGroup` of a :class:`.BuildingBlock`
            placed on that vertex. Only vertices which need to have
            their default alignment changed need to be present in the
            :class:`dict`. If ``None`` then the first :class:`.Edge`
            in :attr:`.Vertex.edges` is used for each vertex. Changing
            which :class:`.Edge` is used will mean that the topology
            graph will represent a different structural isomer.

            The vertices and edges can also be referred to by their
            indices. The vertices are referred to by their index in
            the class attribute :attr:`vertices` while the edges are
            referred to by their index in :attr:`.Vertex.edges`.

        Returns
        -------
        :class:`dict`
            `vertex_alignments` but with indices replaced with the
            vertices and edges they refer to.

        """

        _vertex_alignments = {}
        for v, e in vertex_alignments.items():
            v = self.vertices[v] if isinstance(v, int) else v
            e = v.edges[e] if isinstance(e, int) else e
            _vertex_alignments[v] = e
        return _vertex_alignments

    def _get_instance_vertices(self, vertex_alignments):
        """
        Create the vertices of the topology graph instance.

        Parameters
        ---------
        vertex_alignments : :class:`dict`
            A mapping from a :class:`.Vertex`, in the class
            attribute :attr:`vertices`, to an :class:`.Edge`, in the
            class attribute :attr:`edges`, connected to it.
            The :class:`.Edge` is used to align the first
            :class:`.FunctionalGroup` of a :class:`.BuildingBlock`
            placed on that vertex. Only vertices which need to have
            their default alignment changed need to be present in the
            :class:`dict`. If ``None`` then the first :class:`.Edge`
            in :attr:`.Vertex.edges` is used for each vertex. Changing
            which :class:`.Edge` is used will mean that the topology
            graph will represent a different structural isomer.

        Returns
        -------
        :class:`list`
            A nested :class:`list` which can be indexed as
            ``vertices[x][y][z]``, which will return a :class:`dict`
            for the unit cell at (x, y, z). The :class:`dict` maps
            the vertices in the class attribute :attr:`vertices` to
            the instance clones, for that unit cell.


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
        vertices = it.product(cells, self.vertices)
        for cell, vertex in vertices:
            x, y, z = cell
            clone = vertex.clone(clear_edges=True)
            clone.set_cell(x, y, z)
            clone.aligner_edge = vertex_alignments.get(
                vertex,
                vertex.edges[0]
            )
            # Shift the clone so that it's within the cell.
            shift = 0
            for axis, dim in zip(cell, self._lattice_constants):
                shift += axis * dim
            clone.set_position(clone.get_position()+shift)

            vertex_clones[x][y][z][vertex] = clone
        return vertex_clones

    def _get_instance_edges(self, vertices):
        """
        Create the edges in the topology graph instance.

        Parameters
        ----------
        vertices : :class:`list`
            A nested :class:`list` which can be indexed as
            ``vertices[x][y][z]``, which will return a :class:`dict`
            for the unit cell at (x, y, z). The :class:`dict` maps
            the vertices in the class attribute :attr:`vertices` to
            the instance clones, for that unit cell.

        Returns
        -------
        :class:`tuple` of :class:`.Edge`
            The edges of the topology graph instance.

        """

        edge_clones = []
        # Make a clone for each edge for each unit cell.
        xdim, ydim, zdim = (range(dim) for dim in self._lattice_size)
        cells = it.product(xdim, ydim, zdim)
        edges = it.product(cells, self.edges)
        for cell, edge in edges:
            x, y, z = cell
            # The cell in which edge.vertices[1] is found.
            periodicity = edge.get_periodicity()
            periodic_cell = np.array(cell) + periodicity
            # Wrap around periodic cells, ie those that are less than 0
            # or greater than the lattice size along any dimension.
            dims = zip(periodic_cell, self._lattice_size)
            cell2_x, cell2_y, cell2_z = np.array([
                (dim+max_dim) % max_dim
                for dim, max_dim in dims
            ])
            # Make a vertex map which accounts for the fact that
            # edge.vertices[1] is in cell2.
            v0 = edge.vertices[0]
            v1 = edge.vertices[1]
            vertex_map = {
                v0: vertices[x][y][z][v0],
                v1: vertices[cell2_x][cell2_y][cell2_z][v1]
            }
            # If the edge is not periodic if periodic_cell is did not
            # have to wrap around.
            dims = zip(periodic_cell, self._lattice_size)
            edge_is_not_periodic = all(
                dim >= 0 and dim < max_dim
                for dim, max_dim in dims
            )
            clone = edge.clone(
                vertex_map=vertex_map,
                recalculate_position=edge_is_not_periodic
            )
            edge_clones.append(clone)
            if edge_is_not_periodic:
                clone.set_periodicity(0, 0, 0)
            # Set the aligner edge to the clone.
            for vertex in vertex_map.values():
                if vertex.aligner_edge is edge:
                    vertex.aligner_edge = clone
        return tuple(edge_clones)

    def _before_react(self, mol, vertex_clones, edge_clones):
        if self._periodic:
            return vertex_clones, edge_clones
        return vertex_clones, [
            edge for edge in edge_clones
            if all(dim == 0 for dim in edge.get_periodicity())
        ]

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
            bb = bb_by_degree[len(vertex.edges)]
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
            f'{v.id}: {v.edges.index(v.aligner_edge)}'
            # Only get the vertices in the first unit cell.
            for v in self.vertices[:len(self.__class__.vertices)]
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

    _vertices = (
        _COFVertex(*((1/3)*_a + (1/3)*_b + (1/2)*_c)),
        _COFVertex(*((2/3)*_a + (2/3)*_b + (1/2)*_c))
    )

    vertices = (
        *_vertices,
        _COFVertex.init_at_center(_vertices[0], _vertices[1]),
        _COFVertex.init_at_shifted_center(
            vertices=(_vertices[0], _vertices[1]),
            shifts=((0, 0, 0), (0, -1, 0)),
            lattice_constants=_lattice_constants
        ),
        _COFVertex.init_at_shifted_center(
            vertices=(_vertices[0], _vertices[1]),
            shifts=((0, 0, 0), (-1, 0, 0)),
            lattice_constants=_lattice_constants
        )
    )

    edges = (
        Edge(vertices[2], vertices[0]),
        Edge(vertices[2], vertices[1]),

        Edge(vertices[3], vertices[0]),
        Edge(vertices[3], vertices[1], periodicity=(0, -1, 0)),

        Edge(vertices[4], vertices[0]),
        Edge(vertices[4], vertices[1], periodicity=(-1, 0, 0))
    )


class Hexagonal(COF):
    """
    Represents a hexagonal COF topology graph.

    Building blocks with six and two functional groups are required
    for this topology graph.

    See :class:`.COF` for more details and examples.

    Attributes
    ----------
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

    _vertices = (
        _COFVertex(*((1/4)*_a + (1/4)*_b + (1/2)*_c)),
        _COFVertex(*((1/4)*_a + (3/4)*_b + (1/2)*_c)),
        _COFVertex(*((3/4)*_a + (1/4)*_b + (1/2)*_c)),
        _COFVertex(*((3/4)*_a + (3/4)*_b + (1/2)*_c))
    )

    vertices = (
        *_vertices,
        _COFVertex.init_at_center(_vertices[0], _vertices[1]),
        _COFVertex.init_at_center(_vertices[0], _vertices[2]),
        _COFVertex.init_at_center(_vertices[1], _vertices[2]),
        _COFVertex.init_at_center(_vertices[1], _vertices[3]),
        _COFVertex.init_at_center(_vertices[2], _vertices[3]),
        _COFVertex.init_at_shifted_center(
            vertices=(_vertices[0], _vertices[2]),
            shifts=((0, 0, 0), (-1, 0, 0)),
            lattice_constants=_lattice_constants
        ),
        _COFVertex.init_at_shifted_center(
            vertices=(_vertices[0], _vertices[1]),
            shifts=((0, 0, 0), (0, -1, 0)),
            lattice_constants=_lattice_constants
        ),
        _COFVertex.init_at_shifted_center(
            vertices=(_vertices[0], _vertices[3]),
            shifts=((0, 0, 0), (0, -1, 0)),
            lattice_constants=_lattice_constants
        ),
        _COFVertex.init_at_shifted_center(
            vertices=(_vertices[2], _vertices[1]),
            shifts=((0, 0, 0), (1, -1, 0)),
            lattice_constants=_lattice_constants
        ),
        _COFVertex.init_at_shifted_center(
            vertices=(_vertices[2], _vertices[3]),
            shifts=((0, 0, 0), (0, -1, 0)),
            lattice_constants=_lattice_constants
        ),
        _COFVertex.init_at_shifted_center(
            vertices=(_vertices[1], _vertices[3]),
            shifts=((0, 0, 0), (-1, 0, 0)),
            lattice_constants=_lattice_constants
        ),
        _COFVertex.init_at_shifted_center(
            vertices=(_vertices[3], _vertices[0]),
            shifts=((0, 0, 0), (1, 0, 0)),
            lattice_constants=_lattice_constants
        )
    )

    edges = (
        Edge(vertices[4], vertices[0]),
        Edge(vertices[4], vertices[1]),

        Edge(vertices[5], vertices[0]),
        Edge(vertices[5], vertices[2]),

        Edge(vertices[6], vertices[1]),
        Edge(vertices[6], vertices[2]),

        Edge(vertices[7], vertices[1]),
        Edge(vertices[7], vertices[3]),

        Edge(vertices[8], vertices[2]),
        Edge(vertices[8], vertices[3]),

        Edge(vertices[9], vertices[0]),
        Edge(vertices[9], vertices[2], periodicity=(-1, 0, 0)),

        Edge(vertices[10], vertices[0]),
        Edge(vertices[10], vertices[1], periodicity=(0, -1, 0)),

        Edge(vertices[11], vertices[0]),
        Edge(vertices[11], vertices[3], periodicity=(0, -1, 0)),

        Edge(vertices[12], vertices[2]),
        Edge(vertices[12], vertices[1], periodicity=(1, -1, 0)),

        Edge(vertices[13], vertices[2]),
        Edge(vertices[13], vertices[3], periodicity=(0, -1, 0)),

        Edge(vertices[14], vertices[1]),
        Edge(vertices[14], vertices[3], periodicity=(-1, 0, 0)),

        Edge(vertices[15], vertices[3]),
        Edge(vertices[15], vertices[0], periodicity=(1, 0, 0))
    )


class Square(COF):
    """
    Represents a sqaure COF topology graph.

    Building blocks with four and two functional groups are required
    for this topology graph.

    See :class:`.COF` for more details and examples.

    Attributes
    ----------
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

    _vertices = (
        _COFVertex(*((0.5)*_a + (0.5)*_b + (0.5)*_c)),
    )
    vertices = (
        *_vertices,
        _COFVertex.init_at_shifted_center(
            vertices=(_vertices[0], _vertices[0]),
            shifts=((0, 0, 0), (1, 0, 0)),
            lattice_constants=_lattice_constants
        ),
        _COFVertex.init_at_shifted_center(
            vertices=(_vertices[0], _vertices[0]),
            shifts=((0, 0, 0), (0, 1, 0)),
            lattice_constants=_lattice_constants
        )

    )

    edges = (
        Edge(vertices[1], vertices[0]),
        Edge(vertices[1], vertices[0], periodicity=(1, 0, 0)),
        Edge(vertices[2], vertices[0]),
        Edge(vertices[2], vertices[0], periodicity=(0, 1, 0))
    )


class Kagome(COF):
    """
    Represents a kagome COF topology graph.

    Building blocks with four and two functional groups are required
    for this topology graph.

    See :class:`.COF` for more details and examples.

    Attributes
    ----------
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

    _vertices = (
        _COFVertex(*((1/4)*_a + (3/4)*_b + (0.5)*_c)),
        _COFVertex(*((3/4)*_a + (3/4)*_b + (1/2)*_c)),
        _COFVertex(*((3/4)*_a + (1/4)*_b + (1/2)*_c))
    )

    vertices = (
        *_vertices,
        _COFVertex.init_at_center(_vertices[0], _vertices[1]),
        _COFVertex.init_at_center(_vertices[0], _vertices[2]),
        _COFVertex.init_at_center(_vertices[1], _vertices[2]),
        _COFVertex.init_at_shifted_center(
            vertices=(_vertices[0], _vertices[1]),
            shifts=((0, 0, 0), (-1, 0, 0)),
            lattice_constants=_lattice_constants
        ),
        _COFVertex.init_at_shifted_center(
            vertices=(_vertices[0], _vertices[2]),
            shifts=((0, 0, 0), (-1, 1, 0)),
            lattice_constants=_lattice_constants
        ),
        _COFVertex.init_at_shifted_center(
            vertices=(_vertices[1], _vertices[2]),
            shifts=((0, 0, 0), (0, 1, 0)),
            lattice_constants=_lattice_constants
        )

    )

    edges = (
        Edge(vertices[3], vertices[0]),
        Edge(vertices[3], vertices[1]),

        Edge(vertices[4], vertices[0]),
        Edge(vertices[4], vertices[2]),

        Edge(vertices[5], vertices[1]),
        Edge(vertices[5], vertices[2]),

        Edge(vertices[6], vertices[0]),
        Edge(vertices[6], vertices[1], periodicity=(-1, 0, 0)),

        Edge(vertices[7], vertices[0]),
        Edge(vertices[7], vertices[2], periodicity=(-1, 1, 0)),

        Edge(vertices[8], vertices[1]),
        Edge(vertices[8], vertices[2], periodicity=(0, 1, 0))
    )


class LinkerlessHoneycomb(COF):
    """
    Represents a honeycomb COF topology graph.

    Building blocks with three functional groups are required
    for this topology graph.

    See :class:`.COF` for more details and examples.

    Attributes
    ----------
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

    vertices = (
        _COFVertex(*((1/3)*_a + (1/3)*_b + (1/2)*_c)),
        _COFVertex(*((2/3)*_a + (2/3)*_b + (1/2)*_c))
    )

    edges = (
        Edge(vertices[0], vertices[1]),
        Edge(vertices[0], vertices[1], periodicity=(-1, 0, 0)),
        Edge(vertices[0], vertices[1], periodicity=(0, -1, 0))
    )
