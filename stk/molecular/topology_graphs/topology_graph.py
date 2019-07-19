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


class FGPosition:
    """

    Attributes
    ----------
    func_group : :class:`.FunctionalGroup`
        The functional group placed on the :class:`.VertexPosition`.

    """

    def __init__(self):
        self.func_group = None

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return f'FGPosition()'


class Vertex:
    """
    Represents a vertex of a :class:`.TopologyGraph`.

    Attributes
    ----------
    positions : :class:`tuple` of :class:`.VertexPosition`

    _coord : :class:`numpy.ndarray`
        The position of the vertex.

    """

    def __init__(self, x, y, z, degree):
        self._coord = np.array([x, y, z])
        self.positions = tuple(FGPosition() for i in range(degree))

    @staticmethod
    def _add_fg_position_assignment(fn):

        @wraps(fn)
        def inner(self, building_block):
            r = fn(self, building_block)
            self._assign_fg_positions(building_block)
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
        cls.place_building_block = cls._add_fg_position_assignment(
            cls.place_building_block
        )
        cls.place_building_block = cls._add_position_restoration(
            cls.place_building_block
        )

    def apply_scale(self, scale):
        self._coord *= scale
        return self

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

    def _assign_fg_positions(self, building_block):
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

    Methods
    -------
    :meth:`get_bonded_fgs`

    """

    def __init__(self, *fg_positions):
        self._fg_positions = fg_positions

    def get_bonded_fgs(self):
        """
        Get the functional groups connected by the edge.

        Returns
        -------
        :class:`tuple` of :class:`.FunctionalGroup`
            The functional groups connected by the edge.

        """

        return tuple(p.func_group for p in self._vertex_positions)

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
        vertices = [vertex.clone() for vertex in self._vertices]
        positions = {}
        for clone, vertex in zip(vertices, self._vertices):
            for cp, vp in zip(clone.positions, vertex.positions):
                positions[vp] = cp

        edges = []
        for edge in self._edges:
            positions = (positions[p] for p in edge.positions)
            edges.append(Edge(*positions))

        return vertices, edges

    def _prepare(self, mol):
        mol._conformers.append([])

    def _place_building_blocks(self, mol, vertices):
        bb_map = self._get_bb_map(mol)
        conformer_map = self._get_conformer_map(mol)
        scale = self._get_scale(mol, bb_map, conformer_map)
        for vertex in vertices:
            vertex.apply_scale(scale)

        for vertex in vertices:
            bb = bb_map[vertex]
            conformer_id = conformer_map[vertex]
            coords = vertex.place_building_block(bb, conformer_id)
            mol._conformers[-1].extend(coords)

            if len(mol._conformers) == 1:
                mol.atoms.extend(a.clone() for a in bb.atoms)
                mol.bonds.extend(b.clone() for b in bb.bonds)

    def _get_bonded_fgs(self, mol, edges):
        for edge in edges:
            yield edge.get_bonded_fgs()

    def _get_bb_map(self, mol):
        raise NotImplementedError()

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
