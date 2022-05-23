"""
JSON Serialization Utilities for Molecules
==========================================

"""

from __future__ import annotations

from .....molecular import Atom, Bond

AtomicNumber = int
AtomicCharge = int
AtomJson = tuple[
    AtomicNumber,
    AtomicCharge,
]
AtomId = int
BondOrder = int
BondPeriodicity = tuple[int, int, int]
BondJson = tuple[
    AtomId,
    AtomId,
    BondOrder,
    BondPeriodicity,
]


def atom_to_json(
    atom: Atom,
) -> AtomJson:
    """
    Return a JSON representation of `atom`.

    Parameters:

        atom:
            The atom to serialize.

    Returns:

        A JSON representation of `atom`.

    """

    return (
        atom.get_atomic_number(),
        atom.get_charge(),
    )


def bond_to_json(
    bond: Bond,
) -> BondJson:
    """
    Return a JSON representation of `bond`.

    Parameters:

        bond:
            The bond to serialize.

    Returns:

        A JSON representation of `bond`.

    """

    return (
        bond.get_atom1().get_id(),
        bond.get_atom2().get_id(),
        bond.get_order(),
        bond.get_periodicity(),
    )
