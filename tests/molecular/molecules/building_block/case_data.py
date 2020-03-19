class CaseData:
    """
    A test case.

    Attributes
    ----------
    building_block : :class:`.BuildingBlock`
        The building block to test.

    functional_groups : :class:`tuple` of :class:`.FunctionalGroup`
        The functional groups :attr:`building_block` should be holding.

    core_atoms : :class:`tuple` of :class:`.Atom`
        The core atoms which :attr:`building_block` should be holding.

    placers : :class:`tuple` of :class:`.Atom`
        The *placer* atoms, which :attr:`building_block` should be
        holding.

    """

    def __init__(
        self,
        building_block,
        functional_groups,
        core_atoms,
        placers,
    ):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block to test.

        functional_groups : :class:`tuple` of :class:`.FunctionalGroup`
            The functional groups `building_block` should be holding.

        core_atoms : :class:`tuple` of :class:`.Atom`
            The core atoms which `building_block` should be holding.

        placers : :class:`tuple` of :class:`.Atom`
            The *placer* atoms, which `building_block` should be
            holding.

        """

        self.building_block = building_block
        self.functional_groups = functional_groups
        self.core_atoms = core_atoms
        self.placers = placers
