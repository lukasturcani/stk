"""
Optimizer Utilities
===================

This module defines utilities for optimizers.

"""

from collections import defaultdict

import mchammer as mch


def get_mch_bonds(state):
    """
    Yield the bonds of the :mod:`MCHammer` molecule.

    Parameters
    ----------
    state : :class:`.ConstructionState`
        The state of the molecule under construction.

    Yields
    -------
    :class:`MCHammer.Bond`
        A bond in the molecule.

    """

    for i, bond_infos in enumerate(state.get_bond_infos()):
        ba1 = bond_infos.get_bond().get_atom1().get_id()
        ba2 = bond_infos.get_bond().get_atom2().get_id()
        # Must ensure bond atom id ordering is the same here as in
        # line 30. Therefore, sort here.
        ba1, ba2 = sorted((ba1, ba2))
        yield mch.Bond(id=i, atom_ids=(ba1, ba2))


def get_long_bond_ids(state):
    """
    Yield the ids of the long bonds to optimize.

    Parameters
    ----------
    state : :class:`.ConstructionState`
        The state of the molecule under construction.

    Yields
    -------
    :class:`tuple` of :class:`int`
        A pair of atom ids that identify a bond to be optimized.

    """

    for i, bond_infos in enumerate(state.get_bond_infos()):
        ba1 = bond_infos.get_bond().get_atom1().get_id()
        ba2 = bond_infos.get_bond().get_atom2().get_id()
        # None if for constructed bonds.
        if bond_infos.get_building_block() is None:
            yield sorted((ba1, ba2))


def get_subunits(state):
    """
    Get connected graphs based on building block ids.

    Parameters
    ----------
    state : :class:`.ConstructionState`
        The state of the molecule under construction.

    Returns
    -------
    :class:`.dict`
        The subunits of the molecule split by building block id. Key is
        subunit identifier, Value is :class:`list` of atom ids in
        subunit.

    """

    subunits = defaultdict(list)
    for atom_info in state.get_atom_infos():
        subunits[atom_info.get_building_block_id()].append(
            atom_info.get_atom().get_id()
        )

    return subunits
