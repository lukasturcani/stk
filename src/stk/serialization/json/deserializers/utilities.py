"""
JSON Deserialization Utilities
==============================

"""

from stk.molecular import Atom, AtomInfo, Bond, BondInfo


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
    building_block_atom, = building_block.get_atoms(json[2])

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
            building_blocks[json[0]]
            if json[0] is not None
            else None
        ),
        building_block_id=json[1],
    )
