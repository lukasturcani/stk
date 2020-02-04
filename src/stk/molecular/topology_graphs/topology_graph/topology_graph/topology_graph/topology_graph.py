"""
Adding Topology Graphs
======================

To add a new topology graph a new subclass of :class:`.TopologyGraph`
must be added, which implements its virtual methods. Similarly, new
subclasses of :class:`.VertexData`` and :class:`.Vertex` must also be
made and their virtual methods implemented. When the new subclass of
:class:`.TopologyGraph` is initialized, it must create instances of the
:class:`.VertexData`
subclass, together with :class:`.EdgeData` instances. Once your
topology graph has the vertex and edge data it wants, simply run
:meth:`.TopologyGraph.__init__` and you're done.

When creating a :class:`.VertexData` subclass,
:meth:`~.VertexData.get_vertex` needs to be implemented such that it
returns an instance of your :class:`.Vertex` subclass. If you need to
define a new :meth:`__init__` method for either subclass, you will
also need to implement :meth:`clone` for it.

The :class:`.TopologyGraph` subclass can also create
`construction_stages` if parallel construction needs to be broken down
into separate stages. However,
if this is not the case, then an empty :class:`tuple` can simply be
passed.

Why is both :class:`.VertexData` and :class:`.Vertex` needed?
-------------------------------------------------------------

At first, it may appear that having both :class:`.VertexData` and
:class:`.Vertex` is an unnecessary inconvenience, as when you create
a new :class:`.TopologyGraph` subclass you have to subclass both of
these classes rather than just :class:`.Vertex`. The answer is
related to how these two classes reference other objects in the
:class:`.TopologyGraph`.

:class:`.VertexData` and :class:`.EdgeData` objects keep pointers
to each other in the :attr:`~.VertexData.edges` and
:attr:`~.EdgeData.vertices`. This is extremely convenient for
defining a :class:`.TopologyGraph` because its components can directly
reference each other. However, it poses a significant issue for
serialization. Topology graphs may be highly-cyclic structures
and are therefore they may not possible to serialize with off-the-shelf
serialization tools like :mod:`pickle` or :mod:`dill`. However,
serialization is necessary and fundamental for allowing
parallelization of :class:`.TopologyGraph` construction. The
vertices and edges of the graph have to be serialized and sent to
other cores so that they can place and connect building blocks in
parallel. As a  result, :class:`.VertexData` exists to allow a
convenient definition of a :class:`TopologyGraph`, while
:class:`.Vertex` exists to provide a serializable representation of it.
:class:`.Verex` and :class:`Edge` do not reference other objects
directly, instead they refer to them by their :attr:`.Vertex.id`,
which is used to get an index into :attr:`.TopologyGraph.vertices`
and :attr:`.TopologyGraph.edges`.

:meth:`.VertexData.get_vertex` is used to convert :class:`.VertexData`
into its :class:`.Vertex` counterpart.

"""

from .construction_result import ConstructionResult
from ._construction_state import _ConstructionState


class TopologyGraph:
    """
    Represents topology graphs of :class:`.ConstructedMolecule`.

    The topology graph is an abstract representation of a constructed
    molecule. The vertices indicate where building blocks are placed
    and the edges indicate which building blocks have bonds formed
    between them by the construction process.

    Vertices are responsible for placing the building block molecules.
    By initializing the vertices with different parameters, you can
    alter how they position the building block molecules and therefore
    allow the user to easily specify a different structural isomer.

    Once a building block is placed on a vertex, the functional groups
    on the building block must be mapped to the different edges
    connected to the vertex. The number of functional groups in the
    building block must match the number of edges connected to the
    vertex.

    Once the functional groups are mapped to edges, each edge
    represents a reaction between the functional groups mapped to it.
    Note that more than two functional groups can map to the same edge,
    for cases where you are dealing with something really exotic.
    A :class:`.Reaction` between functional groups is selected based
    on the edges mapped to the edge. A :class:`.Reaction` will
    generally create bonds between the atoms of the functional groups.
    After this you will end up with a :class:`.ConstructedMolecule`.

    """

    def construct(self, vertex_assignments):
        """
        Construct a :class:`.ConstructedMolecule`.

        Parameters
        ----------
        vertex_assignments : :class:`.ConstructedMolecule`
            The :class:`.ConstructedMolecule` instance which needs to
            be constructed.

        Returns
        -------
        :class:`.ConstructionResult`
            The data describing the :class:`.ConstructedMolecule`.

        """

        state = self._get_construction_state(vertex_assignments)
        state = self._before_placement(state)
        state = self._place_building_blocks(state)
        state = self._before_reactions(state)
        state = self._run_reactions(state)
        state = self._clean_up(state)
        return ConstructionResult(state)

    def _get_construction_state(self, vertex_assignments):
        return _ConstructionState()

    def _before_reactions(self, state):
        return state

    def _before_placement(self, state):
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

        return state

    def _place_building_blocks(self, state):
        """
        Place building blocks in `mol` on :attr:`vertices`.

        Parameters
        ----------
        state : :class:`._ConstructionState`
            Holds data necessary to construct the molecule.

        Returns
        -------
        :class:`._ConstructionState`
            The :class:`._ConstructionState` updated to account for
            the placed building blocks.

        """

        raise NotImplementedError()

    def run_reactions(self, state):
        raise NotImplementedError()

    def _clean_up(self, state):
        return state

    def get_adjacency_list(self):
        raise NotImplementedError()

    def __str__(self):
        return repr(self)

    def __repr__(self):
        raise NotImplementedError()
