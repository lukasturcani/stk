import typing

from stk._internal.constructed_molecule import ConstructedMolecule
from stk._internal.topology_graphs.topology_graph.topology_graph import (
    TopologyGraph,
)

T = typing.TypeVar("T", bound=TopologyGraph)


class MoleculeRecord(typing.Generic[T]):
    """
    An abstract base class for molecular records used by the EA.

    Notes:

        You might notice that the public methods of this abstract base
        class are implemented. This is a default implementation provided
        purely for convenience. Subclasses can freely ignore or
        override this implementation.
    """

    def __init__(self, topology_graph: T) -> None:
        """
        Parameters:
            topology_graph:
                The topology graph of a :class:`.ConstructedMolecule`.
        """
        self._molecule = ConstructedMolecule(topology_graph)
        self._topology_graph = topology_graph

    def get_molecule(self) -> ConstructedMolecule:
        """
        Get the molecule held by the record.

        Returns:
            The molecule held by the record.
        """
        return self._molecule

    def get_topology_graph(self) -> T:
        """
        Get the topology graph of the molecule.

        Returns:
            The topology graph.
        """
        return self._topology_graph

    def __repr__(self) -> str:
        return f"MoleculeRecord({self._topology_graph!r})"
