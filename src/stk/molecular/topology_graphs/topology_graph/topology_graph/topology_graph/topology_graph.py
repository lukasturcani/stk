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
:class:`.Vertex` is an unnecssary inconvenience, as when you create
a new :class:`.TopologyGraph` subclass you have to subclass both of
these classes rather than just :class:`.Vertex`. The answer is
related to how these two classes reference other objects in the
:class:`.TopologyGraph`.

:class:`.VertexData` and :class:`.EdgeData` objects keep pointers
to each other in the :attr:`~.VertexData.edges` and
:attr:`~.EdgeData.vertices`. This is extremely convenient for
defining a :class:`.TopologyGraph` because its components can directly
reference each other. However, it poses a significant issue for
serialization. Topology graphs are usually highly-cyclic structures
and are therefore often not possible to serialize with off-the-shelf
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

import numpy as np
import pathos
from collections import namedtuple

from ...utilities import vector_angle






PlacementResult = namedtuple(
    'PlacementResult',
    ['building_block', 'vertex', 'assignments']
)


def _place_building_blocks(vertices, edges):

    def inner(vertex, building_block):
        vertex.place_building_block(building_block, vertices, edges)
        assignments = vertex.assign_func_groups_to_edges(
            building_block=building_block,
            vertices=vertices,
            edges=edges
        )
        return PlacementResult(building_block, vertex, assignments)

    return inner


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

    def __init__(
        self,
        vertex_data,
        edge_data,
        construction_stages,
        num_processes
    ):
        """
        Initialize an instance of :class:`.TopologyGraph`.

        Parameters
        ----------
        vertices : :class:`tuple` of :class:`.VertexData`
            The vertices which make up the graph.

        edges : :class:`tuple` of :class:`.EdgeData`
            The edges which make up the graph.

        construction_stages : :class:`tuple` of :class:`callable`
            A collection of callables, each of which takes a
            :class:`.Vertex` and returns ``True`` or ``False``.
            If the first :class:`callable` is applied to a  vertex in
            `vertices`, that vertex is is part of the first
            construction stage. The second :class:`callable` is then
            applied to all vertices not in the first stage and those
            which return ``True`` belong to the second stage and
            so on.

            Vertices which belong to the same construction stage
            all place building blocks together in parallel, before
            placement is done by any vertices which are part of a later
            stage. This breaks down parallel construction into
            serial stages if synchronization between stages is needed.

            If the topology graph is performing construction serially,
            then all vertices which belong to an earlier stage will
            place their building block before those at a later stage.

        num_processes : :class:`int`
            The number of parallel processes to create during
            :meth:`construct`.

        """

        self._set_data_ids(vertex_data)
        self._set_data_ids(edge_data)
        self.vertices = tuple(
            data.get_vertex() for data in vertex_data
        )
        self.edges = tuple(
            data.get_edge() for data in edge_data
        )
        self._construction_stages = construction_stages
        self._set_stages()
        self._num_processes = num_processes

    def _set_data_ids(self, data):
        for i, data in enumerate(data):
            data.id = i

    def _set_stages(self):
        self._stages = tuple(
            [] for i in range(len(self._construction_stages)+1)
        )
        for vertex in self.vertices:
            placed = False
            for i, stage in enumerate(self._construction_stages):
                if stage(vertex):
                    self._stages[i].append(vertex)
                    placed = True
                    break
            if not placed:
                self._stages[-1].append(vertex)

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


        """

        state = self._get_construction_state(vertex_assignments)
        state = self._prepare(state)
        state = self._place_building_blocks(state)
        state = self._run_reactions(state)
        state = self._clean_up(state)
        return ConstructionResult(state)

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
        :class:`NotImplementedError`
            This is a virtual method which needs to be implemented in
            a subclass.

        """

        raise NotImplementedError()

    def _get_scale(self, mol):
        """
        Get the scale used for vertex and edge positions.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The molecule being constructed.

        Returns
        -------
        :class:`float` or :class:`list` of :class:`float`
            The value by which the position of each :class:`Vertex` and
            is :class:`Edge` is scaled. Can be a single number if all
            axes are scaled by the same amount or a :class:`list` of
            three numbers if each axis is scaled by a different value.

        Raises
        ------
        :class:`NotImplementedError`
            This is a virtual method and needs to be implemented in a
            subclass.

        """

        raise NotImplementedError()

    def _get_vertex_clones(self, mol, scale):
        """
        Yield clones of :attr:`vertices`.

        The order of yielded clones corresponds to the order in
        :attr:`vertices`

        Notes
        -----
        Clones are necessary so that multiple :meth:`construct`
        calls can be done asynchronously and so that the state of the
        original :class:`.Vertex` objects is not
        changed by the construction process.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The molecule being constructed.

        scale : :class:`float` or :class:`list` of :class:`float`
            The value by which the position of each :class:`Vertex` is
            scaled. Can be a single number if all axes are scaled by
            the same amount or a :class:`list` of three numbers if
            each axis is scaled by a different value.

        Yields
        -------
        :class:`.Vertex`
            A vertex clone.

        """

        for vertex in self.vertices:
            yield (
                vertex
                .clone()
                .set_constructed_molecule(mol)
                .apply_scale(scale)
            )

    def _get_edge_clones(self, scale):
        """
        Yield clones of :attr:`edges`.

        The order of yielded edges corresponds to the order in
        :attr:`edges`.

        Parameters
        ----------
        scale : :class:`float` or :class:`list` of :class:`float`
            The value by which the position of each :class:`Edge` is
            scaled. Can be a single number if all axes are scaled by
            the same amount or a :class:`list` of three numbers if
            each axis is scaled by a different value.

        Yields
        -------
        :class:`.Edge`
            An edge clone.

        """

        for edge in self.edges:
            yield edge.clone().apply_scale(scale)

    def _before_react(self, mol, vertices, edges):
        return vertices, edges

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

    def _place_building_blocks(self, mol, vertices, edges):
        """
        Place building blocks in `mol` on :attr:`vertices`.

        Parameters
        ----------
        mol : :class:`.ConstructedMolecule`
            The molecule being constructed.

        vertices : :class:`tuple` of :class:`.Vertex`
            The vertex clones used for construction.

        edges : :class:`tuple` of :class:`.Edge`
            The edge clones used for construction.

        Returns
        -------
        None : :class:`NoneType`

        """

        if self._num_processes == 1:
            return self._place_building_blocks_serial(
                mol=mol,
                vertices=vertices,
                edges=edges
            )
        else:
            return self._place_building_blocks_parallel(
                mol=mol,
                vertices=vertices,
                edges=edges
            )

    def _get_atom_map(self, mol, bb, bb_id):
        atom_map = {}
        for atom in bb.atoms:
            atom_clone = atom.clone()
            atom_clone.id = len(mol.atoms)
            atom_clone.building_block = bb
            atom_clone.building_block_id = bb_id
            atom_map[atom] = atom_clone
            mol.atoms.append(atom_clone)
        return atom_map

    def _assign_func_groups_to_edges(
        self,
        mol,
        bb,
        bb_id,
        edges,
        assignments
    ):
        atom_map = self._get_atom_map(mol, bb, bb_id)
        mol.func_groups.extend(
            fg.clone(atom_map) for fg in bb.func_groups
        )
        num_fgs = len(bb.func_groups)
        for fg_id, edge_id in assignments.items():
            edges[edge_id].assign_func_group(
                func_group=mol.func_groups[-num_fgs+fg_id]
            )
        return atom_map

    def _place_building_blocks_serial(self, mol, vertices, edges):
        bb_id = 0

        vertex_building_blocks = {
            vertex: bb
            for bb, vertices in mol.building_block_vertices.items()
            for vertex in vertices
        }
        # Use a shorter alias.
        counter = mol.building_block_counter
        for stage in self._stages:
            for instance_vertex in stage:
                vertex = vertices[instance_vertex.id]
                bb = vertex_building_blocks[instance_vertex]
                original_coords = bb.get_position_matrix()

                mol._position_matrix.extend(
                    vertex.place_building_block(bb, vertices, edges)
                )
                assignments = vertex.assign_func_groups_to_edges(
                    building_block=bb,
                    vertices=vertices,
                    edges=edges
                )
                atom_map = self._assign_func_groups_to_edges(
                    mol=mol,
                    bb=bb,
                    bb_id=bb_id,
                    edges=edges,
                    assignments=assignments
                )
                # Perform additional, miscellaneous operations.
                vertex.after_assign_func_groups_to_edges(
                    building_block=bb,
                    func_groups=mol.func_groups[-len(bb.func_groups):],
                    vertices=vertices,
                    edges=edges
                )

                bb.set_position_matrix(original_coords)
                mol.bonds.extend(b.clone(atom_map) for b in bb.bonds)
                counter.update([bb])
                bb_id += 1

    def _place_building_blocks_parallel(self, mol, vertices, edges):
        bb_id = 0

        vertex_building_blocks = {
            vertex: bb
            for bb, vertices in mol.building_block_vertices.items()
            for vertex in vertices
        }
        bb_map = {
            bb.get_identity_key(): bb
            for bb in mol.get_building_blocks()
        }
        # Use a shorter alias.
        counter = mol.building_block_counter
        place = _place_building_blocks(vertices, edges)
        with pathos.pools.ProcessPool(self._num_processes) as pool:
            for stage in self._stages:
                verts = []
                bbs = []
                for instance_vertex in stage:
                    verts.append(vertices[instance_vertex.id])
                    bbs.append(vertex_building_blocks[instance_vertex])
                results = pool.map(place, verts, bbs)

                for result in results:
                    result_bb = result.building_block
                    bb = bb_map[result_bb.get_identity_key()]

                    mol._position_matrix.extend(
                        result_bb.get_position_matrix()
                    )
                    atom_map = self._assign_func_groups_to_edges(
                        mol=mol,
                        bb=bb,
                        bb_id=bb_id,
                        edges=edges,
                        assignments=result.assignments
                    )

                    # Perform additional, miscellaneous operations.
                    vertex = vertices[result.vertex.id]
                    num_fgs = len(bb.func_groups)
                    vertex.after_assign_func_groups_to_edges(
                        building_block=result_bb,
                        func_groups=mol.func_groups[-num_fgs:],
                        vertices=vertices,
                        edges=edges
                    )
                    mol.bonds.extend(
                        b.clone(atom_map) for b in bb.bonds
                    )
                    counter.update([bb])
                    bb_id += 1

    def _clean_up(self, mol):
        mol._position_matrix = np.array(mol._position_matrix).T
        for i, atom in enumerate(mol.atoms):
            atom.id = i

    def get_adjacency_list(self):
        raise NotImplementedError()

    def __str__(self):
        return repr(self)

    def __repr__(self):
        raise NotImplementedError()
