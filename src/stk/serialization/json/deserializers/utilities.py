from stk.molecular import Atom, Bond, AtomInfo, BondInfo


def to_atom(id, json):
    return Atom(id, json['atomic_number'], json['charge'])


def to_bond(atoms, json):
    return Bond(
        atom1=atoms[json['atom1']],
        atom2=atoms[json['atom2']],
        order=json['order'],
        periodicity=tuple(json['periodicity']),
    )


def to_atom_info(building_blocks, atom, json):
    return AtomInfo(
        atom=atom,
        building_block=building_blocks[json['building_block']],
        building_block_id=json['building_block_id'],
    )


def to_bond_info(building_blocks, bond, json):
    return BondInfo(
        bond=bond,
        building_block=building_blocks[json['building_block']],
        building_block_id=json['building_block_id'],
    )
