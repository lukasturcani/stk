"""
Defines :class:`.TopologyGraph` and related classes.

.. _`adding topology graphs`:

Extending ``stk``: Adding new topology graphs.
----------------------------------------------

General
.......

A new topology class must be defined. The class must inherit
:class:`Topology`. The new topology class will have to define
the methods :meth:`place_mols` and :meth:`bonded_fgs`. A description
of what these methods should do is given by :meth:`Topology.place_mols`
and :meth:`Topology.bonded_fgs`.

The new class may optionally define the methods :meth:`prepare` and
:meth:`cleanup`. The former performs operations on the molecule
before it is joined up and has atoms deleted via
:meth:`.Reactor.react`. The latter performs any final cleanup
operations on the constructed molecule. For example, converting the end
functional groups of a polymer into hydrogen atoms. See also
:meth:`Topology.cleanup`.

During the construction process, every time a building block is placed
in the :class:`.ConstructedMolecule`, new :class:`FunctionalGroup`
instances must be made, which correspond to the functional groups added
by virtue of adding the building block. These must be added to the
:class:`Reactor` held in :attr:`.Topology.reactor`, specifically into
its :attr:`.Reactor.func_groups` attribute. This means that the
reactor will keep the atom ids in these functional groups up to date
when it deletes atoms. However, note that any functional groups
yielded by :meth:`.Topology.bonded_fgs` are automatically added, so
they do not have to be managed manually. If you do not wish to
automatically add the functional groups into
:attr:`.Reactor.func_groups` you can toggle it in
:attr:`Topology.track_fgs`.

Cages
.....

To add a new cage topology a new class should be created, named
after the topology. This class should inherit :class:`.CageTopology`.
This will give access to various methods which are necessary
for dealing with any cage molecule. See the documenation of
:class:`.CageTopology` for more details.

The new class will only need to have five class attributes added:

    1. a :class:`list` called :attr:`vertices`
    2. a :class:`list` called :attr:`edges`
    3. :attr:`n_windows`, which holds the number of windows the cage
       topology has.
    4. :attr:`n_window_types`, which holds the number of different
       window types. For example, if :attr:`n_window_types` is ``2``,
       then the topology will have two kinds of windows, each with a
       different expected size. Windows of the same type are expected
       to be of the same size.

:attr:`vertices` holds instances of :class:`~.cage.base.Vertex`. Each
instance represents a vertex of a cage and needs to be initialized
with the coordinates of that vertex. Vertices of a cage are where
building blocks of cages are placed.

:attr:`edges` holds instances of the :class:`~.cage.base.Edge`. Each
instance represents an edge of a cage and needs to be initialized
with two instances of :class:`~.cage.base.Vertex`. The
:class:`~.cage.base.Vertex` instances
should be held in :attr:`vertices`, as mentioned above. The two
vertices are the ones which the edge connects. Linkers of cages are
placed on edges. The edge instances automatically derive their
positions from the vertices supplied during initialization.

The vertices need to be positioned such that the center of the
topology is at the origin.


"""

import numpy as np
from functools import wraps
import re

from ..functional_groups import Reactor


class Vertex:
    """
    Represents a vertex of a :class:`.TopologyGraph`.

    Attributes
    ----------
    _coord : :class:`numpy.ndarray`
        The position of the vertex.

    _edges : :class:`list` of :class:`.Edge`
        The edges the :class:`Vertex` is connected to.

    Methods
    -------
    :meth:`apply_scale`
    :meth:`clone`
    :meth:`get_position`
    :meth:`place_building_block`

    """

    def __init__(self, x, y, z):
        self._coord = np.array([x, y, z])
        self._edges = []

    @staticmethod
    def _add_func_group_assignment(fn):
        """
        Add :meth:`_add_func_group_assignment` to `fn`.

        Parameters
        ----------
        fn : :class:`callable`
            A :class:`callable` which is to perform
            :meth:`_add_func_group_assignment` after it is called.
            Usually :meth:`place_building_block`.

        Returns
        -------
        :class:`callable`
            A new :class:`callable` which runs
            :meth:`_add_func_group_assignment` after `fn`.

        """

        @wraps(fn)
        def inner(self, building_block):
            r = fn(self, building_block)
            self._assign_func_groups_to_edges(building_block)
            return r

        return inner

    @staticmethod
    def _add_position_restoration(fn):
        """
        Add building block position restoration after `fn` is called.

        Parameters
        ----------
        fn : :class:`callable`
            A :class:`callable` which is to restore the original
            coordinates of a building block passed to `fn` after
            `fn` is called. Usually :meth:`place_building_block`.

        Returns
        -------
        :class:`callable`
            A new :class:`callable` which restores the original
            coordinates of a building block passed to `fn` after `fn`
            is called.

        """

        @wraps(fn)
        def inner(self, building_block):
            pos_mat = building_block.get_position_matrix()
            r = fn(self, building_block)
            building_block.set_position_matrix(pos_mat)
            return r

        return inner

    def __init_subclass__(cls, **kwargs):
        cls.place_building_block = cls._add_func_group_assignment(
            cls.place_building_block
        )
        cls.place_building_block = cls._add_position_restoration(
            cls.place_building_block
        )

    def apply_scale(self, scale):
        """
        Scale the position of the :class:`.Vertex` by `scale`.

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

    def clone(self):
        """
        Create a clone of the instance.

        Returns
        -------
        :class:`Vertex`
            A clone with the same position and connected to the same
            :class:`.Edge` objects.

        """

        clone = self.__class__(*self._coord)
        clone._edges = list(self._edges)
        return clone

    def get_position(self):
        """
        Return the position of the :class:`Vertex`.

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

    def _assign_func_groups_to_edges(self, building_block):
        """
        Assign functional groups to edges.

        Each :class:`.FunctionalGroup` of the `building_block` needs
        to be assigned to one of the :class:`.Edge` instances in
        :attr:`_edges`.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is needs to have
            functional groups assigned to

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

    def __str__(self):
        return repr(self)

    def __repr__(self):
        x, y, z = self._coord
        cls_name = (
            f'{__package__}.{__name__}.{self.__class__.__name__}'
        )
        # Make sure that the name has all the topology_graph submodule
        # names.
        p = re.compile(r'.*?topology_graphs\.(.*)', re.DOTALL)
        cls_name = p.findall(cls_name)
        return f'{cls_name}({x}, {y}, {z})'


class Edge:
    """
    Represents the edge of a topology graph.

    Attributes
    ----------
    vertices : :class:`tuple` of :class:`.Vertex`
        The vertices which the :class:`Edge` connects.

    _func_groups : :class:`list` of :class:`.FunctionalGroup`
        The functional groups which the edge connects.

    _coord : :class:`numpy.ndarray`
        The position of the edge. It is the centroid the
        :attr:`_vertices`.

    Methods
    -------
    :meth:`assign_func_group`
    :meth:`get_func_groups`

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
        self._func_groups = []

        self._coord = 0
        for i, vertex in enumerate(vertices, 1):
            vertex._edges.append(self)
            self._coord += vertex.get_position()
        self._coord = self._coord / i

    def get_func_groups(self):
        """
        Get the functional groups connected by the edge.

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

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'Edge()'


class TopologyGraph:
    """
    A base class for topology graphs of :class:`.ConstructedMolecule`.

    Attributes
    ----------
    vertices : :class:`.Vertex`
        The vertices which make up the topology graph.

    edges : :class:`.Edge`
        The edges which make up the topology graph.

    processes : :class:`int`
        The number of parallel processes to create during
        :meth:`construct`.

    Methods
    -------
    :meth:`__init__`
    :meth:`construct`

    """

    def __init__(self, vertices, edges, processes=1):
        """
        Initialize an instance of :class:`.TopologyGraph`.

        Parameters
        ----------
        vertices : :class:`.Vertex`
            The vertices which make up the graph.

        edges : :class:`.Edge`
            The edges which make up the graph.

        processes : :class:`int`, optional
            The number of parallel processes to create during
            :meth:`construct`.

        """

        self.vertices = vertices
        self.edges = edges
        self.processes = processes

    def construct(self, mol, building_blocks):
        """
        Construct a :class:`.ConstructedMolecule` conformer.

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

        Note
        ----
        This method will modify
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

        Clones are necessary so that multiple :meth:`construct`
        calls can be done asynchronously and so that the state of the
        original :class:`.Vertex` objects is not
        changed by the construction process.

        Returns
        -------
        :class:`dict`
            A mapping from the original :attr:`vertices` to the clones.

        """

        return {vertex: vertex.clone() for vertex in self.vertices}

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

        if self.processes == 1:
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

        # Use a shorter alias.
        counter = mol.building_block_counter
        for bb, vertices in mol.building_block_vertices.items():
            for vertex in vertices:
                clone = vertex_clones[vertex]
                coords = clone.place_building_block(bb)
                mol._position_matrix.extend(coords)

                atom_map = {}
                for atom in bb.atoms:
                    atom_clone = atom.clone()
                    atom_clone.id = len(mol.atoms)
                    atom_map[atom] = atom_clone
                    atom_clone.building_block = bb
                    atom_clone.building_block_id = counter[bb]
                    mol.atoms.append(atom_clone)

                mol.func_groups.extend(
                    fg.clone(atom_map) for fg in bb.func_groups
                )

                mol.bonds.extend(b.clone() for b in bb.bonds)
                counter.update(bb)

    def _place_building_blocks_parallel(self, mol, vertices):
        raise NotImplementedError()

    def _get_bonded_fgs(self, mol, edges):
        for edge in edges:
            yield edge.get_func_groups()

    def _clean_up(self, mol):
        mol._position_matrix = mol._position_matrix.T
        for i, atom in enumerate(self.atoms):
            atom.id = i

    def __str__(self):
        return repr(self)

    def __repr__(self):
        attrs = ', '.join(
            f'{attr}={val}' for attr, val in vars(self)
            if not attr.startswith('_')
        )
        cls_name = (
            f'{__package__}.{__name__}.{self.__class__.__name__}'
        )
        # Make sure that the name has all the topology_graph submodule
        # names.
        p = re.compile(r'.*?topology_graphs\.(.*)', re.DOTALL)
        cls_name = p.findall(cls_name)
        return f'{cls_name}({attrs})'
