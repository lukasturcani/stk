"""
Constructed Molecule Record
===========================

"""

from .molecule import MoleculeRecord


class ConstructedMoleculeRecord(MoleculeRecord):
    """
    Abstract base class for constructed molecule records of the EA.

    Notes
    -----
    You might notice that the public methods of this abstract base
    class are implemented. This is a default implementation provided
    purely for convenience. Subclasses can freely ignore or
    override this implementation.

    """

    def __init__(self, molecule, topology_graph):
        """
        Initialize a :class:`.ConstructedMoleculeRecord` instance.

        Parameters
        ----------
        molecule : :class:`.ConstructedMolecule`
            The molecule the record holds.

        """

        super().__init__(molecule)
        self._topology_graph = topology_graph

    def get_molecule(self):
        """
        Get the molecule held by the record.

        Returns
        -------
        :class:`.ConstructedMolecule`
            The molecule held by the record.

        """

        return self._molecule

    def get_topology_graph(self):
        """
        Get the topology graph of the constructed molecule.

        Returns
        -------
        :class:`.TopologyGraph`
            The topology graph.

        """

        return self._topology_graph