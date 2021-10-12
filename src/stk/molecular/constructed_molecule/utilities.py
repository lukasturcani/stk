from ..molecular_utilities import get_bond_atom_ids
from ..bond_info import BondInfo


__all__ = (
    'get_bond_info_atom_ids',
)


def get_bond_info_atom_ids(
    bond_info: BondInfo,
) -> list[int]:
    return get_bond_atom_ids(bond_info.get_bond())
