"""
Defines :class:`.TopologyGraph` and related classes.

.. _`adding topology graphs`:

Extending ``stk``: Adding new topology graphs.
----------------------------------------------

To add a new topology graph a new subclass of :class:`.TopologyGraph`
must be added, which implements it's virtual methods. Similarly,
a new subclass of :class:`.Vertex` must also be made and its virtual
methods implemented. When the new subclass of :class:`.TopologyGraph`
is initialized, it must create instances of the :class:`.Vertex`
subclass, together with :class:`.Edge` instances. Once your
topology graph has the vertices and edges it wants, simply run
``super().__init__(vertices, edges, processes)`` and you're done.


"""

import numpy as np
import re
from collections import defaultdict

from ..reactor import Reactor
from ...utilities import vector_theta


class Vertex:
    """
    Represents a vertex in a :class:`.TopologyGraph`.

    Attributes
    ----------
    edges : :class:`list` of :class:`.Edge`
        The edges the :class:`Vertex` is connected to.

    """

    def __init__(self, x, y, z):
        """
        Initialize a :class:`.Vertex`.

        Parameters
        ----------
        x : :class:`float`
            The x coordinate.

        y : :class:`float`
            The y coordinate.

        z : :class:`float`
            The z coordinate.

        """

        self._coord = np.array([x, y, z], dtype=np.dtype('float64'))
        self.edges = []

    def apply_scale(self, scale):
        """
        Scale the position of by `scale`.

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

        self._coord *= scale
        return self

    def clone(self, clear_edges=False):
        """
        Return a clone.

        Parameters
        ----------
        clear_edges : :class:`bool`, optional
            If ``True`` the :attr:`edges` attribute of the clone will
            be empty.

        Returns
        -------
        :class:`Vertex`
            The clone.

        """

        clone = self.__class__.__new__(self.__class__)
        clone._coord = np.array(self._coord)
        clone.edges = [] if clear_edges else list(self.edges)
        return clone

    def get_position(self):
        """
        Return the position.

        Returns
        -------
        :class:`numpy.ndarray`
            The position of the :class:`Vertex`.

        """

        return np.array(self._coord)

    def place_building_block(self, building_block):
        """
        Place `building_block` on the :class:`.Vertex`.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is to be placed on the
            vertex.

        Returns
        -------
        :class:`numpy.nadarray`
            The position matrix of `building_block` after being
            placed.

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method, it needs to be implemented in a
            subclass.

        """

        raise NotImplementedError()

    def assign_func_groups_to_edges(self, building_block, fg_map):
        """
        Assign functional groups to edges.

        Each :class:`.FunctionalGroup` of the `building_block` needs
        to be associated with one of the :class:`.Edge` instances in
        :attr:`edges`. Then, using `fg_map`, the
        :class:`FunctionalGroup` instances in the molecule being
        constructed need to be assigned to those edges. This is
        because bonds need to be formed between functional groups of
        the molecule being constructed, not the `building_block`.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is needs to have
            functional groups assigned to

        fg_map : :class:`dict`
            A mapping from :class:`.FunctionalGroup` instances in
            `building_block` to the equivalent
            :class:`.FunctionalGroup` instances in the molecule being
            constructed.

        Returns
        -------
        None : :class:`NoneType`

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method, it needs to be implemented in a
            subclass.

        """

        raise NotImplementedError()

    def _get_edge_centroid(self, edge_ids=None):
        """
        Return the centroid of the connected edges.

        Parameters
        ----------
        edge_ids : :class:`iterable` of :class:`int`
            The ids of edges which are used to calculate the centroid.
            If ``None``, then all  the edges in :attr:`edges` are used.

        Returns
        -------
        :class:`numpy.ndarray`
            The centroid of the edges.

        """

        if edge_ids is None:
            edge_ids = range(len(self.edges))

        edge_positions = []
        for i, edge_id in enumerate(edge_ids, 1):
            edge_positions.append(self.edges[edge_id].get_position())
        return np.sum(edge_positions, axis=0) / i

    def _get_edge_plane_normal(self, edge_ids=None):
        """
        Get the normal to the plane on which the :attr:`edges` lie.

        Parameters
        ----------
        edge_ids : :class:`iterable` of :class:`int`
            The ids of edges which are used to calculate the plane.
            If there are more than three, a plane of best fit across
            edges is returned. If ``None``, then all  the edges in
            :attr:`edges` are used.

        Returns
        -------
        :class:`numpy.ndarray`
            A unit vector which describes the normal to the plane of
            the edges.

        Raises
        ------
        :class:`ValueError`
            If there are not at least 3 edges, which is necessary to
            define a plane.

        """

        if edge_ids is None:
            edge_ids = range(len(self.edges))
        else:
            # The iterable is used mutliple times.
            edge_ids = list(edge_ids)

        if len(edge_ids) < 3:
            raise ValueError(
                'At least 3 edges are necessary to create a plane.'
            )

        edge_positions = []
        for i, edge_id in enumerate(edge_ids, 1):
            edge_positions.append(self.edges[edge_id].get_position())
        edge_positions = np.array(edge_positions)

        centroid = np.sum(edge_positions, axis=0) / i
        normal = np.linalg.svd(edge_positions - centroid)[-1][2, :]
        if vector_theta(normal, centroid) > np.pi/2:
            normal *= -1
        return normal

    def __str__(self):
        return repr(self)

    def __repr__(self):
        x, y, z = self._coord
        cls_name = (
            f'{__name__}.{self.__class__.__name__}'
        )
        # Make sure that the name has all the topology_graph submodule
        # names.
        p = re.compile(r'.*?topology_graphs\.(.*)', re.DOTALL)
        cls_name = p.findall(cls_name)[0]
        return f'{cls_name}({x}, {y}, {z})'


class Edge:
    """
    Represents an edge in a topology graph.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which the :class:`Edge` connects.

    """

    def __init__(self, *vertices):
        """
        Initialize an :class:`Edge`.

        Parameters
        ----------
        *vertices : :class:`.Vertex`
            The vertices which the :class:`Edge` connects.

        """

        self.vertices = vertices
        # The FunctionalGroup instances which the edge connects.
        # These will belong to the molecules placed on the vertices
        # connected by the edge.
        self._func_groups = []

        self._coord = 0
        for i, vertex in enumerate(vertices, 1):
            vertex.edges.append(self)
            self._coord += vertex.get_position()
        self._coord = self._coord / i

    def get_func_groups(self):
        """
        Get the functional groups connected by this edge.

        Returns
        -------
        :class:`tuple` of :class:`.FunctionalGroup`
            The functional groups connected by the edge.

        """

        return tuple(self._func_groups)

    def assign_func_group(self, func_group):
        """
        Assign `func_group` to be connected by this edge.

        Parameters
        ----------
        func_group : :class:`.FunctionalGroup`
            The functional group to be assigned to the edge.

        Returns
        -------
        :class:`Edge`
            The edge is returned.

        """

        self._func_groups.append(func_group)

    def get_position(self):
        """
        Return the position.

        Returns
        -------
        :class:`numpy.ndarray`
            The position of the :class:`Edge`.

        """

        return np.array(self._coord)

    def __str__(self):
        return repr(self)

    def __repr__(self):
        vertices = ', '.join(repr(v) for v in self.vertices)
        cls_name = (
            f'{__name__}.{self.__class__.__name__}'
        )
        # Make sure that the name has all the topology_graph submodule
        # names.
        p = re.compile(r'.*?topology_graphs\.(.*)', re.DOTALL)
        cls_name = p.findall(cls_name)[0]
        return f'{cls_name}({vertices})'


class TopologyGraph:
    """
    Represents topology graphs of :class:`.ConstructedMolecule`.

    The topology graph is an abstract representation of a constructed
    molecule. The vertices indicate where building blocks are placed
    and the edges indicate which building blocks have bonds formed
    between them by the construction process.

    Vertices are responsible for placing the building block molecules.
    By initializing the vertices with different settings, they can
    position the building block molecules differently and therefore
    allow the user to easily specify a different structural isomer.

    Once a building block is placed on a vertex, the functional groups
    on the building block must be assigned to the different edges
    connected to the vertex. The number of functional groups in the
    building block must match the number of edges connected to the
    vertex.

    Once the functional groups are assigned to edges, each edge
    represents a reaction between the functional groups assigned to it.
    Note that an edge can be assigned more than two functional groups,
    in case you are dealing with something really exotic. The
    functional groups are then matched to an appropriate reaction,
    which generally creates bonds between the atoms of the functional
    groups. After this you will end up with a
    :class:`.ConstructedMolecule`.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`tuple` of :class:`.Edge`
        The edges which make up the topology graph.

    """

    def __init__(self, vertices, edges, processes):
        """
        Initialize an instance of :class:`.TopologyGraph`.

        Parameters
        ----------
        vertices : :class:`tuple` of :class:`.Vertex`
            The vertices which make up the graph.

        edges : :class:`tuple` of :class:`.Edge`
            The edges which make up the graph.

        processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        """

        self.vertices = vertices
        self.edges = edges
        self._processes = processes

    def construct(self, mol, building_blocks):
        """
        Construct a :class:`.ConstructedMolecule`.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The :class:`.ConstructedMolecule` instance which needs to
            be constructed.

        building_blocks : :class:`list` of :class:`.Molecule`
            The :class:`.BuildingBlock` and
            :class:`ConstructedMolecule` instances which
            represent the building block molecules used for
            construction. Only one instance is present per building
            block molecule, even if multiples of that building block
            join up to form the :class:`ConstructedMolecule`.

        Returns
        -------
        None : :class:`NoneType`

        """

        if mol.building_block_vertices is None:
            mol.building_block_vertices = defaultdict(list)
            self._assign_building_blocks_to_vertices(
                mol=mol,
                building_blocks=building_blocks
            )

        vertex_clones = self._clone_vertices()
        edge_clones = self._clone_edges(vertex_clones)

        self._prepare(mol)
        self._place_building_blocks(mol, vertex_clones)

        reactor = Reactor(mol)
        for fgs in self._get_bonded_fgs(mol, edge_clones):
            reactor.add_reaction(*fgs)
        mol.bonds_made = reactor.finalize()

        self._clean_up(mol)

    def _assign_building_blocks_to_vertices(
        self,
        mol,
        building_blocks
    ):
        """
        Assign `building_blocks` to :attr:`vertices`.

        Assignment is done by modifying
        :attr:`.ConstructedMolecule.building_block_vertices`.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The :class:`.ConstructedMolecule` instance being
            constructed.

        building_blocks : :class:`list` of :class:`.Molecule`
            The :class:`.BuildingBlock` and
            :class:`ConstructedMolecule` instances which
            represent the building block molecules used for
            construction. Only one instance is present per building
            block molecule, even if multiples of that building block
            join up to form the :class:`ConstructedMolecule`.

        Returns
        -------
        None : :class:`NoneType`

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method which needs to be implemented in
            a subclass.

        """

        raise NotImplementedError()

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

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method and needs to be implemented in a
            subclass.

        """

        raise NotImplementedError()

    def _clone_vertices(self):
        """
        Create clones of :attr:`vertices`.

        Notes
        -----
        Clones are necessary so that multiple :meth:`construct`
        calls can be done asynchronously and so that the state of the
        original :class:`.Vertex` objects is not
        changed by the construction process.

        Returns
        -------
        :class:`dict`
            A mapping from the original :attr:`vertices` to the clones.

        """

        return {
            vertex: vertex.clone(clear_edges=True)
            for vertex in self.vertices
        }

    def _clone_edges(self, vertex_clones):
        """
        Create clones of :attr:`edges`.

        Parameters
        ----------
        vertex_clones : :class:`dict`
            A mapping from the original :attr:`vertices` to the
            clones.

        Returns
        -------
        :class:`list` of :class:`.Edge`
            The cloned :attr:`edges`.

        """

        edges = []
        for edge in self.edges:
            vertices = (
                vertex_clones[vertex] for vertex in edge.vertices
            )
            edges.append(Edge(*vertices))
        return edges

    def _prepare(self, mol):
        """
        Do preprocessing on `mol` before construction.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The molecule being constructed.

        Returns
        -------
        None : :class:`NoneType`

        """

        return

    def _place_building_blocks(self, mol, vertex_clones):
        """
        Place building blocks in `mol` on :attr:`vertices`.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The molecule being constructed.

        vertex_clones : :class:`dict`
            A mapping from :attr:`vertices` to their clones being used
            for a particular :meth:`construct` call.

        Returns
        -------
        None : :class:`NoneType`

        """

        if self._processes == 1:
            return self._place_building_blocks_serial(
                mol=mol,
                vertex_clones=vertex_clones
            )
        else:
            return self._place_building_blocks_parallel(
                mol=mol,
                vertex_clones=vertex_clones
            )

    def _place_building_blocks_serial(self, mol, vertex_clones):
        scale = self._get_scale(mol)
        for vertex in vertex_clones.values():
            vertex.apply_scale(scale)

        bb_id = 0
        # Use a shorter alias.
        counter = mol.building_block_counter
        for bb, vertices in mol.building_block_vertices.items():
            original_coords = bb.get_position_matrix()

            for vertex in vertices:
                # Use a clone of the vertex for positioning so that
                # state of the originals is not changed.
                vertex_clone = vertex_clones[vertex]
                # Get the coordinates of the building block when
                # placed on the vertex clone.
                coords = vertex_clone.place_building_block(bb)
                # Add the coordinates to the constructed molecule.
                mol._position_matrix.extend(coords)

                # Create a map from each atom in the building block to
                # a clone of that atom in the constructed molecule.
                atom_map = {}
                for atom in bb.atoms:
                    atom_clone = atom.clone()
                    atom_clone.id = len(mol.atoms)
                    atom_map[atom] = atom_clone

                    atom_clone.building_block = bb
                    atom_clone.building_block_id = bb_id

                    mol.atoms.append(atom_clone)

                # Create a map from each functional group in the
                # building block to a clone of that functional group
                # in the constructed molecule.
                fg_map = {
                    fg: fg.clone(atom_map) for fg in bb.func_groups
                }
                mol.func_groups.extend(fg_map.values())

                # Assign the functional groups in the contructed
                # molecule to edges in the topology graph.
                vertex_clone.assign_func_groups_to_edges(bb, fg_map)

                bb.set_position_matrix(original_coords)
                mol.bonds.extend(b.clone(atom_map) for b in bb.bonds)
                counter.update([bb])
                bb_id += 1

    def _place_building_blocks_parallel(self, mol, vertices):
        raise NotImplementedError('TODO')

    def _get_bonded_fgs(self, mol, edges):
        for edge in edges:
            yield edge.get_func_groups()

    def _clean_up(self, mol):
        mol._position_matrix = np.array(mol._position_matrix).T
        for i, atom in enumerate(mol.atoms):
            atom.id = i

    def __str__(self):
        return repr(self)
