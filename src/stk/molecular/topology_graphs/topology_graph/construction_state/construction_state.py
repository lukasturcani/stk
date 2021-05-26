"""
Construction State
==================

"""

from .graph_state import _GraphState
from .molecule_state import _MoleculeState


class ConstructionState:
    """
    The state of the molecule and topology graph under construction.

    """

    def __init__(
        self,
        building_block_vertices,
        edges,
        lattice_constants=(),
    ):
        """
        Initialize a :class:`.ConstructionState` instance.

        Parameters
        ----------
        building_block_vertices : :class:`dict`
            Maps each :class:`.BuildingBlock` to be placed, to a
            :class:`tuple` of :class:`.Vertex` instances, on which
            it should be placed.

        edges : :class:`tuple` of :class:`.Edge`
            The edges of the topology graph.

        lattice_constants : :class:`tuple`, optional
            A :class:`numpy.ndarray` for each lattice constant.
            Can be an empty :class:`tuple` if the topology graph is
            not periodic.

        """

        self._graph_state = _GraphState(
            building_block_vertices=building_block_vertices,
            edges=edges,
            lattice_constants=lattice_constants,
        )
        self._molecule_state = _MoleculeState()

    def clone(self):
        """
        Return a clone.

        Returns
        -------
        :class:`.ConstructionState`
            The clone. Has the same type as the original instance.

        """

        clone = self.__class__.__new__(self.__class__)
        clone._graph_state = self._graph_state
        clone._molecule_state = self._molecule_state
        return clone

    def _with_placement_results(
        self,
        vertices,
        edges,
        building_blocks,
        results,
    ):
        """
        Modify the instance.

        """

        self._molecule_state = (
            self._molecule_state.with_placement_results(
                vertices=vertices,
                edges=edges,
                building_blocks=building_blocks,
                results=results,
            )
        )
        return self

    def with_placement_results(
        self,
        vertices,
        edges,
        building_blocks,
        results,
    ):
        """
        Return a clone holding the placement results.

        Parameters
        ----------
        vertices : :class:`tuple` of :class:`.Vertex`
            The vertices used for placement.

        edges : :class:`tuple`
            For each vertex in `vertices`, a :class:`tuple` of
            :class:`.Edge` instances connected to it.

        building_blocks : :class:`tuple` of :class:`.BuildingBlock`
            For each vertex in `vertices`, the building block placed
            on it.

        results : :class:`tuple` of :class:`._PlacementResult`
            For every vertex in `vertices`, the result of the
            placement.

        Returns
        -------
        :class:`.ConstructionState`
            The clone holding the placement results. Has the same
            type as the original instance.

        """

        return self.clone()._with_placement_results(
            vertices=vertices,
            edges=edges,
            building_blocks=building_blocks,
            results=results,
        )

    def get_lattice_constants(self):
        """
        Get the lattice constants of the state.

        Returns
        -------
        :class:`tuple` of :class:`numpy.ndarray`
            The lattice constants.

        """

        return self._graph_state.get_lattice_constants()

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

        return self._graph_state.get_building_block(vertex_id)

    def get_vertices(self, vertex_ids):
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

        yield from self._graph_state.get_vertices(vertex_ids)

    def get_num_vertices(self):
        """
        Get the number of vertices in the topology graph.

        Returns
        -------
        :class:`int`
            The number of vertices in the topology graph.

        """

        return self._graph_state.get_num_vertices()

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

        return self._graph_state.get_edge(edge_id)

    def get_num_edges(self):
        """
        Get the number of edges in the topology graph.

        Returns
        -------
        :class:`int`
            The number of edges.

        """

        return self._graph_state.get_num_edges()

    def get_edges(self, vertex_id):
        """
        Get the edges connect to a vertex.

        Parameters
        ----------
        vertex_id : :class:`int`
            The id of a vertex.

        Returns
        -------
        :class:`tuple` of :class:`.Edge`
            The connected edges.

        """

        return self._graph_state.get_edges(vertex_id)

    def get_edge_group_functional_groups(self, edge_group):
        """
        Yield the functional groups associated with `edge_group`.

        Parameters
        ----------
        edge_group : :class:`.EdgeGroup`
            The edge group, whose functional groups are desired.

        Yields
        ------
        :class:`.FunctionalGroup`
            A functional group which belongs to `edge_group`.

        """

        yield from (
            self._molecule_state.get_edge_group_functional_groups(
                edge_group=edge_group,
            )
        )

    def _with_reaction_results(self, reactions, results):
        """
        Modify the instance.

        """

        self._molecule_state = (
            self._molecule_state.with_reaction_results(
                reactions=reactions,
                results=results,
            )
        )
        return self

    def with_reaction_results(self, reactions, results):
        """
        Return a clone holding the reaction results.

        Parameters
        ----------
        reactions : :class:`tuple` of :class:`.Reaction`
            The reactions.

        results : :class:`.ReactionResult`
            For each reaction in `reactions`, its result.

        Returns
        -------
        :class:`.ConstructionState`
            The clone holding the reaction results. Has the same type
            as the original instance.

        """

        return self.clone()._with_reaction_results(reactions, results)

    def _with_lattice_constants(self, lattice_constants):
        """
        Modify the instance.

        """

        self._graph_state = (
            self._graph_state.with_lattice_constants(
                lattice_constants=lattice_constants,
            )
        )
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
        :class:`.ConstructionState`
            The clone holding the new lattice constants. Has the same
            type as the original instance.

        """

        return self.clone()._with_lattice_constants(lattice_constants)

    def _with_position_matrix(self, position_matrix):
        """
        Modify the instance.

        """

        self._molecule_state = (
            self._molecule_state.with_position_matrix(
                position_matrix=position_matrix,
            )
        )
        return self

    def with_position_matrix(self, position_matrix):
        """
        Return a clone holding the `position_matrix`.

        Parameters
        ----------
        position_matrix : :class:`numpy.ndarray`
            The position matrix of the clone. The shape of the matrix
            is ``(n, 3)``.

        Returns
        -------
        :class:`.ConstructionState`
            The clone holding the new position matrix. Has the same
            type as the original instance.

        """

        return self.clone()._with_position_matrix(position_matrix)

    def _with_vertices(self, vertices):
        """
        Modify the instance.

        """

        self._graph_state = self._graph_state.with_vertices(vertices)
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
        :class:`.ConstructionState`
            The clone. Has the same type as the original instance.

        """

        return self.clone()._with_vertices(vertices)

    def get_position_matrix(self):
        """
        Get the position matrix of the molecule being constructed.

        Returns
        -------
        :class:`numpy.ndarray`
            The position matrix.

        """

        return self._molecule_state.get_position_matrix()

    def get_atoms(self):
        """
        Yield the atoms of the molecule being constructed.

        Yields
        ------
        :class:`.Atom`
            An atom of the molecule being constructed.

        """

        yield from self._molecule_state.get_atoms()

    def get_bonds(self):
        """
        Yield the bonds of the molecule being constructed.

        Yields
        ------
        :class:`.Bond`
            A bond of the molecule being constructed.

        """

        yield from self._molecule_state.get_bonds()

    def get_atom_infos(self):
        """
        Yield the atom infos of the molecule being constructed.

        Yields
        ------
        :class:`.AtomInfo`
            An atom info of the molecule being constructed.

        """

        yield from self._molecule_state.get_atom_infos()

    def get_bond_infos(self):
        """
        Yield the bond infos of the molecule being constructed.

        Yields
        ------
        :class:`.BondInfo`
            The bond info of the molecule being constructed.

        """

        yield from self._molecule_state.get_bond_infos()

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

        return self._graph_state.get_num_building_block(building_block)

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

        yield from self._graph_state.get_building_blocks()
