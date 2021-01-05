"""
This module defines utilities for optimizers.

"""

import mchammer as mch
from collections import defaultdict


class OptimizationIncompleteError(Exception):
    ...


def get_mch_bond_topology(state):
    """
    Returns bonds with atom1_id < atom2_id and long bonds to optimize.

    """

    long_bond_ids = []
    mch_bonds = []
    for i, bond_infos in enumerate(state.get_bond_infos()):
        ba1 = bond_infos.get_bond().get_atom1().get_id()
        ba2 = bond_infos.get_bond().get_atom2().get_id()
        # Must ensure bond atom id ordering is such that
        # atom1_id < atom2_id.
        if ba1 < ba2:
            mch_bonds.append(
                Bond(id=i, atom1_id=ba1, atom2_id=ba2)
            )
        else:
            mch_bonds.append(
                Bond(id=i, atom1_id=ba2, atom2_id=ba1)
            )

        # Define long bonds based on bond_info.
        # None for constructed bonds.
        if bond_infos.get_building_block() is None:
            if ba1 < ba2:
                ids = (ba1, ba2)
            else:
                ids = (ba2, ba1)
            long_bond_ids.append(ids)

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
