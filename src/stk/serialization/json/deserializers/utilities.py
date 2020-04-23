from stk.molecular import Atom, Bond, AtomInfo, BondInfo


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
    return AtomInfo(
        atom=atom,
        building_block=(
            building_blocks[json[0]]
            if json[0] is not None
            else None
        ),
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
