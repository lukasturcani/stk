class CaseData:
    """
    A test case.

    Attributes
    ----------
    vertex : :class:`.Vertex`
        The vertex to test.

    edges : :class:`tuple` of :class:`.Edge`
        The edges connected to :attr:`.vertex`

    building_block : :class:`.BuildingBlock`
        The building block to be placed by :attr:`.vertex`.

    position : :class:`numpy.ndarray`
        The correct centroid of atoms in :attr:`.position_ids` after
        placement.

    alignment_tests : :class:`dict`
        Maps a :class:`callable` to a :class:`numpy.ndarray`.
        Each :class:`callable` takes a singe parameter,
        `building_block`, and returns a :class:`numpy.ndarray`.

        Each :class:`callable` in `alignment_tests` represents a
        test, and the array which it maps to is the correct result
        for that test. If the array returned by the :class:`callable`
        does not match the array it is mapped to, the test will fail.

        The point of these tests is to make sure that
        :attr:`.building_block` is aligned correctly after placement.
        Therefore, alignment tests should return some vector which
        depends on the specific orientation of the building block.

    functional_group_edges : :class:`dict`
        The correct mapping of functional group id to edge id.

    position_ids : :class:`tuple` of :class:`int`
        The ids of atoms which are used to determine if the building
        block was positioned correctly.

    """

    def __init__(
        self,
        vertex,
        edges,
        building_block,
        position,
        alignment_tests,
        functional_group_edges,
        position_ids=None,
    ):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        vertex : :class:`.Vertex`
            The vertex to test.

        edges : :class:`tuple` of :class:`.Edge`
            The edges connected to `vertex`

        building_block : :class:`.BuildingBlock`
            The building block to be placed by `vertex`.

        position : :class:`numpy.ndarray`
            The correct centroid of atoms in `position_ids` after
            placement.

        alignment_tests : :class:`dict`
            Maps a :class:`callable` to a :class:`numpy.ndarray`.
            Each :class:`callable` takes a singe parameter,
            `building_block`, and returns a :class:`numpy.ndarray`.

            Each :class:`callable` in `alignment_tests` represents a
            test, and the array which it maps to is the correct result
            for that test. If the array returned by the
            :class:`callable` does not match the array it is mapped to,
            the test will fail.

            The point of these tests is to make sure that
            `.building_block` is aligned correctly after
            placement. Therefore, alignment tests should return some
            vector which depends on the specific orientation of the
            building block.

        functional_group_edges : :class:`dict`
            The correct mapping of functional group id to edge id.

        position_ids : :class:`tuple` of :class:`int`, optional
            The ids of atoms which are used to determine if the
            building block was positioned correctly. If ``None``, the
            placer ids will be used.

        """

        if position_ids is None:
            position_ids = tuple(building_block.get_placer_ids())

        self.vertex = vertex
        self.edges = edges
        self.building_block = building_block
        self.position = position
        self.alignment_tests = alignment_tests
        self.functional_group_edges = functional_group_edges
        self.position_ids = position_ids

    def __str__(self):
        return (
            f"CaseData({self.vertex}, {self.edges}, "
            f"{self.building_block}, {self.position}, "
            f"{self.alignment_tests}, {self.functional_group_edges})"
        )

    def __repr__(self):
        return str(self)
