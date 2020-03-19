import stk


class CaseData:
    """
    A test case.

    Attributes
    ----------
    bond_info : :class:`.BondInfo`
        The bond info to test.

    bond : :class:`.Bond`
        The bond, which should be held by :attr:`.bond_info`.

    building_block : :class:`.BuildingBlock`
        The building block, which should be held by :attr:`.bond_info`.

    building_block_id : :class:`int`
        The building block id, which should be held by
        :attr:`.bond_info`.

    """

    def __init__(self, bond, building_block, building_block_id):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        bond : :class:`.Bond`
            The bond, which should be held by a bond info.

        building_block : :class:`.BuildingBlock`
            The building block, which should be held by a bond info.

        building_block_id : :class:`int`
            The building block id, which should be held by a bond info.

        """

        self.bond_info = stk.BondInfo(
            bond=bond,
            building_block=building_block,
            building_block_id=building_block_id,
        )
        self.bond = bond
        self.building_block = building_block
        self.building_block_id = building_block_id
