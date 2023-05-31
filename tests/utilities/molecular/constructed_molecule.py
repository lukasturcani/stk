import itertools as it

from .molecule import (
    is_equivalent_atom,
    is_equivalent_bond,
    is_equivalent_molecule,
)


def is_equivalent_constructed_molecule(
    constructed_molecule1,
    constructed_molecule2,
):
    is_equivalent_molecule(
        molecule1=constructed_molecule1,
        molecule2=constructed_molecule2,
    )
    are_equivalent_atom_infos(
        infos1=constructed_molecule1.get_atom_infos(),
        infos2=constructed_molecule2.get_atom_infos(),
    )
    are_equivalent_bond_infos(
        infos1=constructed_molecule1.get_bond_infos(),
        infos2=constructed_molecule2.get_bond_infos(),
    )

    for building_block1, building_block2 in it.zip_longest(
        constructed_molecule1.get_building_blocks(),
        constructed_molecule2.get_building_blocks(),
    ):
        is_equivalent_molecule(
            building_block1,
            building_block2,
        )

    for building_block1, building_block2 in zip(
        constructed_molecule1.get_building_blocks(),
        constructed_molecule2.get_building_blocks(),
    ):
        num1 = constructed_molecule1.get_num_building_block(
            building_block=building_block1,
        )
        num2 = constructed_molecule2.get_num_building_block(
            building_block=building_block2,
        )
        assert num1 == num2


def are_equivalent_atom_infos(infos1, infos2):
    for info1, info2 in it.zip_longest(infos1, infos2):
        is_equivalent_atom_info(info1, info2)


def are_equivalent_bond_infos(infos1, infos2):
    for info1, info2 in it.zip_longest(infos1, infos2):
        is_equivalent_bond_info(info1, info2)


def is_equivalent_atom_info(info1, info2):
    is_equivalent_atom(info1.get_atom(), info2.get_atom())
    if info1.get_building_block() is None:
        assert info2.get_building_block() is None
    else:
        is_equivalent_molecule(
            info1.get_building_block(),
            info2.get_building_block(),
        )

    if info1.get_building_block_atom() is None:
        assert info2.get_building_block_atom() is None
    else:
        is_equivalent_atom(
            atom1=info1.get_building_block_atom(),
            atom2=info2.get_building_block_atom(),
        )

    assert info1.get_building_block_id() == info2.get_building_block_id()


def is_equivalent_bond_info(info1, info2):
    is_equivalent_bond(info1.get_bond(), info2.get_bond())
    if info1.get_building_block() is None:
        assert info2.get_building_block() is None
    else:
        is_equivalent_molecule(
            info1.get_building_block(),
            info2.get_building_block(),
        )
    assert info1.get_building_block_id() == info2.get_building_block_id()
