import stk


class CaseData:
    """
    A test case.

    Attributes
    ----------
    atom_info : :class:`.AtomInfo`
        The atom info to test.

    atom : :class:`.Atom`
        The atom, which should be held by :attr:`.atom_info`.

    building_block_atom : :class:`.Atom`
        The building block atom which should be held by
        :attr:`.atom_info`.

    building_block : :class:`.BuildingBlock`
        The building block, which should be held by :attr:`.atom_info`.

    building_block_id : :class:`int`
        The building block id, which should be held by
        :attr:`.atom_info`.

    """

    def __init__(
        self,
        atom,
        building_block_atom,
        building_block,
        building_block_id,
    ):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        atom : :class:`.Atom`
            The atom, which should be held by an :class:`.AtomInfo`.

        building_block_atom : :class:`.Atom`
            The building block atom which should be held by an
            :class:`.AtomInfo`.

        building_block : :class:`.BuildingBlock`
            The building block, which should be held by an
            :class:`.AtomInfo`.

        building_block_id : :class:`int`
            The building block id, which should be held by an
            :class:`.AtomInfo`.

        """

        self.atom_info = stk.AtomInfo(
            atom=atom,
            building_block_atom=building_block_atom,
            building_block=building_block,
            building_block_id=building_block_id,
        )
        self.atom = atom
        self.building_block_atom = building_block_atom
        self.building_block = building_block
        self.building_block_id = building_block_id
