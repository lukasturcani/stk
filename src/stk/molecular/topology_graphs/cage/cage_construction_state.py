from ..topology_graph import ConstructionState


class _CageConstructionState(ConstructionState):
    def __init__(
        self,
        building_block_vertices,
        edges,
        num_placement_stages,
        vertex_degrees,
        lattice_constants,
    ):
        super().__init__(
            building_block_vertices=building_block_vertices,
            edges=edges,
            lattice_constants=lattice_constants,
        )
        self._num_placement_stages = num_placement_stages
        self._num_placement_stages_done = 0
        self._vertex_degrees = dict(vertex_degrees)
        self._neighbor_positions = {}

    def _with_placement_results(
        self,
        vertices,
        edges,
        building_blocks,
        results,
    ):
        # Need to iterate multiple times through results.
        results = tuple(results)
        super()._with_placement_results(
            vertices=vertices,
            edges=edges,
            building_blocks=building_blocks,
            results=results,
        )

        # No need to update vertex positions after the last stage.
        self._num_placement_stages_done += 1
        if (
            self._num_placement_stages_done
            >= self._num_placement_stages
        ):
            return self

        self._update_neighbor_positions(
            vertices=vertices,
            edges=edges,
            building_blocks=building_blocks,
            results=results,
        )
        return self

    def _update_neighbor_positions(
        self,
        vertices,
        edges,
        building_blocks,
        results,
    ):
        for vertex, vertex_edges, building_block, result in zip(
            vertices,
            edges,
            building_blocks,
            results,
        ):
            building_block = building_block.with_position_matrix(
                position_matrix=result.position_matrix,
            )
            edge_functional_groups = dict(zip(
                result.functional_group_edges.values(),
                result.functional_group_edges.keys(),
            ))
            for neighbor_id, edge_id in self._get_neighbors(
                vertex=vertex,
                vertex_edges=vertex_edges,
            ):
                fg_id = edge_functional_groups[edge_id]
                functional_group = next(
                    building_block.get_functional_groups(fg_id)
                )
                self._neighbor_positions[neighbor_id] = (
                    self._neighbor_positions.get(neighbor_id, [])
                )
                self._neighbor_positions[neighbor_id].append(
                    building_block.get_centroid(
                        atom_ids=functional_group.get_placer_ids(),
                    )
                )

        self._graph_state = self._graph_state.with_vertices(
            vertices=self._get_new_vertices(),
        )

    def _get_neighbors(self, vertex, vertex_edges):
        for edge in vertex_edges:
            neighbor_id = (
                edge.get_vertex1_id()
                if vertex.get_id() != edge.get_vertex1_id()
                else edge.get_vertex2_id()
            )
            neighbor, = self._graph_state.get_vertices(neighbor_id)
            if neighbor.use_neighbor_placement():
                yield neighbor_id, edge.get_id()

    def _get_new_vertices(self):
        for vertex_id, vertex in enumerate(
            self._graph_state.get_vertices()
        ):
            if (
                len(self._neighbor_positions.get(vertex_id, []))
                == self._vertex_degrees[vertex_id]
            ):
                yield vertex.with_position(
                    position=(
                        sum(self._neighbor_positions[vertex_id])
                        / len(self._neighbor_positions[vertex_id])
                    ),
                )
            else:
                yield vertex

    def clone(self):
        clone = super().clone()
        clone._neighbor_positions = {
            key: list(value)
            for key, value in self._neighbor_positions.items()
        }
        clone._num_placement_stages_done = (
            self._num_placement_stages_done
        )
        clone._num_placement_stages = self._num_placement_stages
        clone._vertex_degrees = dict(self._vertex_degrees)
        return clone
