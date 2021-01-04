"""
This module defines utilities for optimizers.

"""

from mchammer import Bond

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


def merge_subunits(state, subunits):
    """
    Merge subunits in stk.Molecule by building block ids.

    """

    subunit_building_block_ids = {i: set() for i in subunits}
    atom_infos = list(state.get_atom_infos())
    for su in subunits:
        su_ids = subunits[su]
        for su_id in su_ids:
            atom_info = [
                i for i in atom_infos
                if i.get_atom().get_id() == su_id
            ][0]

            subunit_building_block_ids[su].add(
                atom_info.get_building_block_id()
            )

    new_subunits = {}
    taken_subunits = set()
    for su in subunits:
        bb_ids = subunit_building_block_ids[su]
        if len(bb_ids) > 1:
            raise ValueError(
                'Subunits not made up of singular BuildingBlock'
            )
        bb_id = list(bb_ids)[0]
        if su in taken_subunits:
            continue

        compound_subunit = subunits[su]
        has_same_bb_id = [
            (su_id, bb_id) for su_id in subunits
            if list(subunit_building_block_ids[su_id])[0] == bb_id
            and su_id != su
        ]

        for su_id, bb_id in has_same_bb_id:
            for i in subunits[su_id]:
                compound_subunit.add(i)
            taken_subunits.add(su_id)
        new_subunits[su] = compound_subunit

    return new_subunits
