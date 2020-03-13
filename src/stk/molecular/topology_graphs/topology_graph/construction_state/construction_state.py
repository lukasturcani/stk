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

        edges : :class:`tuple` of :class:`.Edge`
            The edges of the topology graph.

        lattice_constants : :class:`tuple`, optional
            A :class:`numpy.ndarray` for each lattice constants.
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
        return self.clone()._with_placement_results(
            vertices=vertices,
            edges=edges,
            building_blocks=building_blocks,
            results=results,
        )

    def get_building_block(self, vertex_id):
        """

        """

        return self._graph_state.get_building_block(vertex_id)

    def get_vertex(self, vertex_id):
        """

        """

        return self._graph_state.get_vertex(vertex_id)

    def get_num_vertices(self):
        return self._graph_state.get_num_vertices()

    def get_edge(self, edge_id):
        return self._graph_state.get_edge(edge_id)

    def get_num_edges(self):
        return self._graph_state.get_num_edges()

    def get_edges(self, vertex_id):
        return self._graph_state.get_edges(vertex_id)

    def get_edge_group_functional_groups(self, edge_group):
        yield from (
            self._molecule_state.get_edge_group_functional_groups(
                edge_group=edge_group,
            )
        )

    def _with_reaction_results(self, reactions, results):
        self._molecule_state = (
            self._molecule_state.with_reaction_results(
                reactions=reactions,
                results=results,
            )
        )
        return self

    def with_reaction_results(self, reactions, results):
        return self.clone()._with_reaction_results(reactions, results)

    def _with_vertices(self, vertices):
        self._graph_state = self._graph_state.with_vertices(vertices)
        return self

    def with_vertices(self, vertices):
        return self.clone()._with_vertices(vertices)

    def get_position_matrix(self):
        """

        """

        return self._molecule_state.get_position_matrix()

    def get_atoms(self):
        yield from self._molecule_state.get_atoms()

    def get_bonds(self):
        yield from self._molecule_state.get_bonds()

    def get_atom_infos(self):
        yield from self._molecule_state.get_atom_infos()

    def get_bond_infos(self):
        yield from self._molecule_state.get_bond_infos()
