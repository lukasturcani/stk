import typing

from stk._internal.constructed_molecule import ConstructedMolecule
from stk._internal.topology_graphs.topology_graph.topology_graph import (
    TopologyGraph,
)


class MoleculeRecord:
    """
    An abstract base class for molecular records used by the EA.

    Notes:

        You might notice that the public methods of this abstract base
        class are implemented. This is a default implementation provided
        purely for convenience. Subclasses can freely ignore or
        override this implementation.

    """

    def __init__(self, topology_graph: TopologyGraph) -> None:
        """
        Parameters:
            topology_graph:
                The topology graph of a :class:`.ConstructedMolecule`.
        """
        self._molecule = ConstructedMolecule(topology_graph)
        self._topology_graph = topology_graph
        self._fitness_value = None
        self._normalized_fitness_value = None

    def get_molecule(self) -> ConstructedMolecule:
        """
        Get the molecule held by the record.

        Returns:
            The molecule held by the record.

        """

        return self._molecule

    def get_topology_graph(self) -> TopologyGraph:
        """
        Get the topology graph of the molecule.

        Returns:
            The topology graph.

        """

        return self._topology_graph

    def get_fitness_value(
        self, normalized: bool = True
    ) -> typing.Any | float | None:
        """
        Get the fitness value of the molecule in the record.

        Parameters:

            normalized:
                Toggles the return of the normalized vs unnormalized
                fitness value. The unnormalized fitness value is
                guaranteed to be constant for the same molecule
                across generations, while the normalized one is allowed
                to change.

        Returns:
            The fitness value.
        """

        return (
            self._normalized_fitness_value
            if normalized
            else self._fitness_value
        )

    def clone(self) -> typing.Self:
        """
        Return a clone.

        Returns:
            MoleculeRecord: The clone. Has the same type as
            the original record.

        """

        clone = self.__class__.__new__(self.__class__)
        clone._molecule = self._molecule
        clone._topology_graph = self._topology_graph
        clone._fitness_value = self._fitness_value
        clone._normalized_fitness_value = self._normalized_fitness_value
        return clone

    def with_fitness_value(
        self,
        fitness_value: typing.Any,
        normalized: bool = True,
    ) -> typing.Self:
        """
        Return a clone holding a different fitness value.

        Parameters:
            fitness_value:
                The fitness value of the clone.

            normalized:
                Toggles if the normalized or unnormalized fitness value is
                being set. If ``False``, both the normalized and
                unnormalized fitness values with be set to `fitness_value`.

        Returns:
            MoleculeRecord: The clone. Has the same type as
            the original record.

        """

        return self.clone()._with_fitness_value(
            fitness_value=fitness_value,
            normalized=normalized,
        )

    def _with_fitness_value(self, fitness_value, normalized):
        if not normalized:
            self._fitness_value = fitness_value
        self._normalized_fitness_value = fitness_value
        return self
