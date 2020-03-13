from ......bonds import BondInfo


class _BondBatchData:
    """
    Holds data about a batch of bonds.

    """

    __slots__ = ['_bonds', '_bond_infos']

    def __init__(self, bonds, atom_map):
        """
        Initialize a :class:`.BondBatchData` instance.

        Parameters
        ----------
        bonds : :class:`iterable` of :class:`.Bond`
            The bonds, for which data should be created.

        atom_map : :class:`dict`
            Maps the ids of atoms held by `bonds`, to the new atoms,
            which the bonds in the :class:`.BondBatchData` instance
            should hold.

        """

        self._bonds = _bonds = []
        self._bond_infos = bond_infos = []

        for bond in bonds:
            _bonds.append(bond.with_atoms(atom_map))
            bond_infos.append(BondInfo(_bonds[-1], None, None))

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
