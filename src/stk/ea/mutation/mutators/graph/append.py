import numpy as np

from .mutator import GraphMutator
from ...records import GraphMutationRecord
from ....molecule_records import MolecularGraphRecord


class GraphAppend(GraphMutator):
    """

    """

    def __init__(
        self,
        building_blocks,
        name='GraphAppend',
        random_seed=None,
    ):
        self._building_blocks = building_blocks
        self._name = name
        self._generator = np.random.RandomState(random_seed)

    def mutate(self, record):
        graph = record.get_topology_graph()
        free_nodes = tuple(
            node for node in graph.get_nodes()
            if node.get_num_functional_groups() != 0
        )
        # Need some error handling here in case len(free_nodes) is 0.
        extended_node = self._generator.choice(free_nodes)
        new_graph = Graph(
            nodes=_append(
                iterable=map(_to_node, graph.get_nodes()),
                last=Node(
                    building_block=...,
                ),
            ),
            edges=_append(
                iterable=map(_to_edge, graph.get_edges()),
                last=Edge(
                    functional_group1=FunctionalGroupReference(
                        node_id=extended_node.get_id(),
                        functional_group_id=...,
                    ),
                    functional_group2=FunctionalGroupReference(
                        node_id=graph.get_num_nodes(),
                        functional_group_id=...,
                    ),
                ),
            ),
            reaction_factory=...,
        )

        return GraphMutationRecord(
            molecule_record=MolecularGraphRecord(new_graph),
            mutator_name=self._name,
        )


def _append(iterable, last):
    yield from iterable
    yield last


def _to_node(node):
    return Node(
        building_block=node.get_building_block(),
        position=node.get_position(),
    )
