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
    :meth:`get_position`
    :meth:`clone`

    """

    def __init__(self, x, y, z, degree):
        self._coord = np.array([x, y, z])
        self._edges = []

    @staticmethod
    def _add_func_group_assignment(fn):

        @wraps(fn)
        def inner(self, building_block):
            r = fn(self, building_block)
            self._assign_func_groups_to_edges(building_block)
            return r

        return inner

    @staticmethod
    def _add_position_restoration(fn):

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
        self._coord *= scale
        return self

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
        Place a building block molecule on the :class:`.Vertex`.

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is to be placed on the
            vertex.

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method, it needs to be implemented in a
            subclass.

        """

        raise NotImplementedError()

    def _assign_func_groups_to_edges(self, building_block):
        """
        Assign

        Parameters
        ----------
        building_block : :class:`.Molecule`
            The building block molecule which is needs to have
            functional groups assigned to

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
        return f'Vertex({x}, {y}, {z}, degree={self.degree})'


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

    Methods
    -------
    :meth:`__init__`
    :meth:`construct`

    """

    def __init__(self, vertices, edges):
        """
        Initialize an instance of :class:`.TopologyGraph`.

        Parameters
        ----------
        vertices : :class:`.Vertex`
            The vertices which make up the graph.

        edges : :class:`.Edge`
            The edges which make up the graph.

        """

        self.vertices = vertices
        self.edges = edges

    def construct(self, mol):
        """
        Construct a :class:`.ConstructedMolecule` conformer.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The :class:`.ConstructedMolecule` instance which needs to
            be constructed.

        Returns
        -------
        None : :class:`NoneType`

        """

        vertices, edges = self._clone_vertices_and_edges()

        self._reactor = Reactor()
        self._place_building_blocks(mol, vertices)
        self._prepare(mol)

        self._reactor.set_molecule(mol.mol)
        mol.func_groups = self._reactor.func_groups

        for fgs in self._get_bonded_fgs(mol, edges):
            self._reactor.react(*fgs, track_fgs=self._track_fgs)
        mol.mol = self._reactor.result(self._del_atoms)
        mol.bonds_made = self._reactor.bonds_made

        self._clean_up(mol)

        # Reactor can't be pickled because it contains an EditableMol,
        # which can't be pickled.
        self._reactor = None

    def _get_scale(self, mol, bb_map, conformer_map):
        raise NotImplementedError()

    def _clone_vertices_and_edges(self, mol):
        vertex_clones = {
            vertex: vertex.clone() for vertex in self.vertices
        }
        edges = []
        for edge in self.edges:
            vertices = (
                vertex_clones[vertex] for vertex in edge.vertices
            )
            edges.append(Edge(*vertices))
        return list(vertex_clones.values()), edges

    def _prepare(self, mol):
        return

    def _place_building_blocks(self, mol, vertices):
        scale = self._get_scale(mol, mol.bb_map)
        for vertex in vertices:
            vertex.apply_scale(scale)

        for vertex in vertices:
            bb = mol.bb_map[vertex]
            coords = vertex.place_building_block(bb)
            mol._conformers[-1].extend(coords)

            if len(mol._conformers) == 1:
                mol.atoms.extend(a.clone() for a in bb.atoms)
                mol.bonds.extend(b.clone() for b in bb.bonds)

    def _get_bonded_fgs(self, mol, edges):
        for edge in edges:
            yield edge.get_func_groups()

    def _clean_up(self, mol):
        mol._conformers[-1] = mol.conformers[-1].T

    def __str__(self):
        return repr(self)

    def __repr__(self):
        attrs = ', '.join(
            f'{attr}={val}' for attr, val in vars(self)
            if not attr.startswith('_')
        )
        return f'{self.__class__.__name__}({attrs})'
