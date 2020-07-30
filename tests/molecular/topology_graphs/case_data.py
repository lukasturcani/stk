class CaseData:
    """
    A :class:`.TopologyGraph` test case.

    Attributes
    ----------
    topology_graph : :class:`.TopologyGraph`
        The topology graph to test.

    cell : :class:`tuple` of :class:`np.array` or :class:`NoneType`
        Tuple of cell lattice vectors (shape: (3,)) in Angstrom.
        `None` if topology graph is not periodic.

    """

    def __init__(self, topology_graph, cell):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        topology_graph : :class:`.TopologyGraph`
            The topology graph to test.

        cell : :class:`tuple` of :class:`np.array` or :class:`NoneType`
            Tuple of cell lattice vectors (shape: (3,)) in Angstrom.
            `None` if topology graph is not periodic.

        """

        self.topology_graph = topology_graph
        self.cell = cell
