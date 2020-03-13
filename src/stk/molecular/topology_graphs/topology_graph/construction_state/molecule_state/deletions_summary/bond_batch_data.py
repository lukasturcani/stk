from .....bonds import BondInfo


class _BondBatchData:
    __slots__ = ['_bonds', '_bond_infos']

    def __init__(
        self,
        bonds,
        atom_map,
        building_block=None,
        building_block_id=None,
    ):
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
        yield from self._bonds

    def get_bond_infos(self):
        yield from self._bond_infos
