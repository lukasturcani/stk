"""
This module defines utilities for optimizers.

"""

import mchammer as mch
from collections import defaultdict


def get_mch_bond_topology(state):
    """
    Returns bond topology including long bonds to optimize.

    """

    long_bond_ids = []
    mch_bonds = []
    for i, bond_infos in enumerate(state.get_bond_infos()):
        ba1 = bond_infos.get_bond().get_atom1().get_id()
        ba2 = bond_infos.get_bond().get_atom2().get_id()
        # Must ensure bond atom id ordering is the same here as in
        # line 30. Therefore, sort here.
        ba1, ba2 = sorted((ba1, ba2))
        mch_bonds.append(mch.Bond(id=i, atom_ids=(ba1, ba2)))

        # Define long bonds based on bond_info.
        # None for constructed bonds.
        if bond_infos.get_building_block() is None:
            long_bond_ids.append((ba1, ba2))

    return long_bond_ids, mch_bonds


def get_subunits(state):
    """
    Get connected graphs based on building block ids.

    Returns
    -------
    subunits : :class:`.dict`
        The subunits of `mol` split by building block id. Key is
        subunit identifier, Value is :class:`iterable` of atom ids in
        subunit.

    """

    subunits = defaultdict(list)
    for atom_info in state.get_atom_infos():
        subunits[atom_info.get_building_block_id()].append(
            atom_info.get_atom().get_id()
        )

    return subunits
