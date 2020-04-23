from ......bonds import BondInfo


class _BondBatch:
    """
    A batch of bonds.

    """

    __slots__ = ['_bonds', '_bond_infos']

    def __init__(
        self,
        bonds,
        atom_map,
        building_block,
        building_block_id,
    ):
        """
        Initialize a :class:`.BondBatch` instance.

        Parameters
        ----------
        bonds : :class:`iterable` of :class:`.Bond`
            The bonds, which should be added to the batch.

        atom_map : :class:`dict`
            Maps the ids of atoms held by `bonds`, to the new atoms,
            which the bonds in the :class:`.BondBatchData` instance
            should hold.

        building_block : :class:`.BuildingBlock`
            The building block from which the bonds originate.

        building_block_id : :class:`.int`
            An id, unique to that building block and placement.

        """

        self._bonds = _bonds = []
        self._bond_infos = bond_infos = []

        for bond in bonds:
            _bonds.append(bond.with_atoms(atom_map))
            bond_infos.append(
                BondInfo(
                    bond=_bonds[-1],
                    building_block=building_block,
                    building_block_id=building_block_id,
                )
            )

    def get_bonds(self):
        """
        Yield the bonds in the batch.

        Yields
        ------
        :class:`.Bond`
            A bond in the batch.

        """

        yield from self._bonds

    def get_bond_infos(self):
        """
        Yield info about the bonds in the batch.

        Yields
        ------
        :class:`.BondInfo`
            Info about a bond in the batch.

        """

        yield from self._bond_infos
