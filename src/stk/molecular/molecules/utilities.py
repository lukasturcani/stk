"""
Molecular Utilities
===================

This module defines utilities for molecules.

"""


def sort_bond_atoms_by_id(bond):
    if bond.get_atom1().get_id() < bond.get_atom2().get_id():
        return bond
    elif bond.get_order() == 9:
        # Bond ordering in Dative bonds cannot change.
        return bond
    else:
        return bond.__class__(
            atom1=bond.get_atom2(),
            atom2=bond.get_atom1(),
            order=bond.get_order(),
            periodicity=bond.get_periodicity(),
        )


def get_bond_atom_ids(bond):
    return sorted((
        bond.get_atom1().get_id(),
        bond.get_atom2().get_id(),
    ))


def get_bond_info_atom_ids(bond_info):
    return get_bond_atom_ids(bond_info.get_bond())
