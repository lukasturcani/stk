"""
Defines the base :class:`Topology` type.

.. _`adding topologies`:

Extending stk: Adding new topologies.
-------------------------------------

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

from inspect import signature
from collections import Counter
import numpy as np
from functools import wraps

from ..functional_groups import Reactor


class TopologyMeta(type):
    """
    Makes a repr of an instance, based initialization arguments used.

    """

    def __call__(self, *args, **kwargs):

        # Get the arguments, keyword arguments and defulat initialized
        # arguments used to make an instance of Topology.
        sig = signature(self.__init__).bind_partial(self,
                                                    *args, **kwargs)
        sig.apply_defaults()
        sig = dict(sig.arguments)
        sig.pop('self')
        # Create the Topology instance.
        obj = super().__call__(*args, **kwargs)
        # Use the arguments the object was initialized with to make
        # a repr of the object and place it in the `repr` attribute.
        # The __repr__() function in Topology will then just return
        # this attribute.
        c = ', '.join(
            f'{key!s}={value!r}' for key, value in sorted(sig.items())
        )
        obj._repr = f'{self.__name__}({c})'
        return obj


class Topology(metaclass=TopologyMeta):
    """
    Constructs :class:`.ConstructedMolecule` from building blocks.

    More accurately, child classes of :class:`Topology` take care of
    constructing :class:`.ConstructedMolecule` from building blocks.

    This class directly defines any operations and attributes that are
    needed by any topology during construction. However, this class is
    not used directly. It is intended to be inherited. All
    :attr:`.ConstructedMolecule.topology` attributes hold an instance
    of a :class:`Topology` child class. Child classes of
    :class:`Topology` define operations specific to that one topology.
    For example, each child class must define a :meth:`join_mols`,
    which creates bonds between the building blocks of a
    :class:`.ConstructedMolecule`. The way in which this is done will
    depend on what kind of molecules are being constructed. In
    addition, each child class must define methods which place the
    building blocks in appropriate positions.

    Attributes
    ----------
    _del_atoms : :class:`bool`
        Toggles whether deleter atoms are deleted by
        :meth:`.Reactor.result`.

    _reactor : :class:`.Reactor`
        The reactor which performs the reactions.

    _track_fgs : :class:`bool`
        Toggles whether functional groups yielded by
        :meth:`bonded_fgs` are automatically added into
        :attr:`.Reactor.func_groups`.

    Methods
    -------
    :meth:`construct`

    """

    def __init__(self, del_atoms=True, track_fgs=True):
        self._del_atoms = del_atoms
        self._track_fgs = track_fgs

    def construct(self, mol):
        raise NotImplementedError()

    def __str__(self):
        return repr(self)

    def __repr__(self):
        # The `_repr` attribute is made in the TopologyMeta __call__()
        # method, when the Topology object is instantiated.
        return self._repr

    def __eq__(self, other):
        return repr(self) == repr(other)

    def __hash__(self):
        return id(self)


class VertexPosition:
    def __init__(self):
        self.func_group = None


class Vertex:
    """

    """

    def __init__(self, x, y, z, degree):
        self._coord = np.array([x, y, z])
        self.positions = [VertexPosition() for i in range(degree)]

    @staticmethod
    def _add_vertex_position_assignment(fn):

        @wraps(fn)
        def inner(self, bb, conformer_id):
            r = fn(self, bb, conformer_id)
            self._assign_vertex_positions(bb, conformer_id)
            return r

        return inner

    @staticmethod
    def _add_position_restoration(fn):

        @wraps(fn)
        def inner(self, bb, conformer_id):
            pos_mat = bb.get_position_matrix(conformer_id=conformer_id)
            r = fn(self, bb, conformer_id)
            bb.set_position_matrix(pos_mat, conformer_id)
            return r

        return inner

    def __init_subclass__(cls, **kwargs):
        cls.place_building_block = cls._add_vertex_position_assignment(
            cls.place_building_block
        )
        cls.place_building_block = cls._add_position_restoration(
            cls.place_building_block
        )

    def apply_scale(self, scale):
        self._coord *= scale
        return self

    def place_building_block(self, bb, conformer_id):
        raise NotImplementedError()

    def _assign_vertex_positions(self, bb, conformer_id):
        raise NotImplementedError()


class Edge:
    def __init__(self, *vertex_positions):
        self._vertex_positions = vertex_positions

    def get_bonded_fgs(self):
        return [p.func_group for p in self._vertex_positions]


class TopologyGraph(Topology):
    """

    """

    def __init__(self, vertices, edges):
        self._vertices = vertices
        self._edges = edges

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

        mol.bonds_made = 0
        mol.bb_counter = Counter()

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

    def _get_conformer_map(self, mol):
        raise NotImplementedError()

    def _clean_up(self, mol):
        mol._conformers[-1] = mol.conformers[-1].T
