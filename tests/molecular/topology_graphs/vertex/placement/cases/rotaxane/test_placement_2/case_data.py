class CaseData:
    """
    A test case.

    Attributes
    ----------
    vertex1 : :class:`.Vertex`
        The first vertex to test.

    vertex2 : :class:`.Vertex`
        The second vertex to test.

    building_block : :class:`.BuildingBlock`
        The building block placed by the vertices.

    atom_ids : :class:`tuple` of :class:`int`
        The ids of atoms used for calculating the building block
        orientation.

    """

    def __init__(self, vertex1, vertex2, building_block, atom_ids):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        vertex1 : :class:`.Vertex`
            The first vertex to test.

        vertex2 : :class:`.Vertex`
            The second vertex to test.

        building_block : :class:`.BuildingBlock`
            The building block placed by the vertices.

        atom_ids : :class:`tuple` of :class:`int`
            The ids of atoms used for calculating the building block
            orientation.

        """

        self.vertex1 = vertex1
        self.vertex2 = vertex2
        self.building_block = building_block
        self.atom_ids = atom_ids
