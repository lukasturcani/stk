"""
Molecular Utilities
===================

This module defines utilities for molecules.

"""

from .bond import Bond


__all__ = (
    'sort_bond_atoms_by_id',
    'get_bond_atom_ids',
)


def sort_bond_atoms_by_id(bond: Bond) -> Bond:
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


def get_bond_atom_ids(bond: Bond) -> list[int]:
    return sorted((
        bond.get_atom1().get_id(),
        bond.get_atom2().get_id(),
    ))
