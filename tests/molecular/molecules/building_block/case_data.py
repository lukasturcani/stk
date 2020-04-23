class CaseData:
    """
    A test case.

    Attributes
    ----------
    building_block : :class:`.BuildingBlock`
        The building block to test.

    functional_groups : :class:`tuple` of :class:`.FunctionalGroup`
        The functional groups :attr:`.building_block` should be
        holding.

    core_atom_ids : :class:`tuple` of :class:`int`
        The correct core atom ids for :attr:`.building_block`.

    placer_ids : :class:`tuple` of :class:`int`
        The correct *placer* ids for :attr:`.building_block`.

    """

    def __init__(
        self,
        building_block,
        functional_groups,
        core_atom_ids,
        placer_ids,
    ):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block to test.

        functional_groups : :class:`tuple` of :class:`.FunctionalGroup`
            The functional groups `building_block` should be holding.

        core_atom_ids : :class:`tuple` of :class:`int`
            The correct core atom ids for `building_block`.

        placer_ids : :class:`tuple` of :class:`int`
            The correct *placer* ids for `building_block`.

        """

        self.building_block = building_block
        self.functional_groups = functional_groups
        self.core_atom_ids = core_atom_ids
        self.placer_ids = placer_ids
