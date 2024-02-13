import stk


class CaseData:
    """
    A test case.

    Attributes:
        building_block:
            The building block, which will be written to a file, so
            that it can be initialized from it.

        functional_groups:
            The functional groups :attr:`.building_block` should be
            holding.

        known_repr:
            The representation of the building block.

        core_atom_ids:
            The correct core atom ids for :attr:`.building_block`.

        placer_ids:
            The correct *placer* ids for :attr:`.building_block`.

    """

    def __init__(
        self,
        building_block: stk.BuildingBlock,
        functional_groups: tuple,
        known_repr: str,
        core_atom_ids: tuple[int, ...],
        placer_ids: tuple[int, ...],
    ) -> None:
        """
        Initialize a :class:`.CaseData` instance.

        Parameters:
            building_block:
                The building block, which will be written to a file, so
                that it can be initialized from it.

            functional_groups:
                The functional groups :attr:`.building_block` should be
                holding.

            known_repr:
                The representation of the building block.

            core_atom_ids:
                The correct core atom ids for :attr:`.building_block`.

            placer_ids:
                The correct *placer* ids for :attr:`.building_block`.

        """

        self.building_block = building_block
        self.functional_groups = functional_groups
        self.known_repr = known_repr
        self.core_atom_ids = core_atom_ids
        self.placer_ids = placer_ids
