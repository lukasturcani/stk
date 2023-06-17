"""
JSON Deserialization Utilities
==============================

"""

from stk._internal.atom import Atom
from stk._internal.atom_info import AtomInfo
from stk._internal.bond import Bond
from stk._internal.bond_info import BondInfo


def to_atom(id, json):
    return Atom(id, json[0], json[1])


def to_bond(atoms, json):
    return Bond(
        atom1=atoms[json[0]],
        atom2=atoms[json[1]],
        order=json[2],
        periodicity=tuple(json[3]),
    )


def to_atom_info(building_blocks, atom, json):
    if json[0] is None:
        return AtomInfo(
            atom=atom,
            building_block_atom=None,
            building_block=None,
            building_block_id=None,
        )

    building_block = building_blocks[json[0]]
    (building_block_atom,) = building_block.get_atoms(json[2])

    return AtomInfo(
        atom=atom,
        building_block_atom=building_block_atom,
        building_block=building_block,
        building_block_id=json[1],
    )


def to_bond_info(building_blocks, bond, json):
    return BondInfo(
        bond=bond,
        building_block=(
            building_blocks[json[0]] if json[0] is not None else None
        ),
        building_block_id=json[1],
    )


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
