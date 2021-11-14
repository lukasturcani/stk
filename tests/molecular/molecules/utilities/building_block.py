import itertools as it

from tests.utilities import is_equivalent_atom, is_equivalent_molecule


def are_equivalent_functional_groups(
    functional_groups1,
    functional_groups2,
):
    functional_groups = it.zip_longest(
        functional_groups1,
        functional_groups2,
    )
    for fg1, fg2 in functional_groups:
        is_equivalent_functional_group(fg1, fg2)


def is_equivalent_functional_group(
    functional_group1,
    functional_group2,
):
    assert functional_group1.__class__ is functional_group2.__class__

    atoms = it.zip_longest(
        functional_group1.get_atoms(),
        functional_group2.get_atoms(),
    )
    for atom1, atom2 in atoms:
        is_equivalent_atom(atom1, atom2)

    for placer_id1, placer_id2 in it.zip_longest(
        functional_group1.get_placer_ids(),
        functional_group2.get_placer_ids(),
    ):
        assert placer_id1 == placer_id2

    for core_atom_id1, core_atom_id2 in it.zip_longest(
        functional_group1.get_core_atom_ids(),
        functional_group2.get_core_atom_ids(),
    ):
        assert core_atom_id1 == core_atom_id2


def is_equivalent_building_block(building_block1, building_block2):
    is_equivalent_molecule(building_block1, building_block2)
    are_equivalent_functional_groups(
        functional_groups1=building_block1.get_functional_groups(),
        functional_groups2=building_block2.get_functional_groups(),
    )
    for placer_id1, placer_id2 in it.zip_longest(
        building_block1.get_placer_ids(),
        building_block2.get_placer_ids(),
    ):
        assert placer_id1 == placer_id2

    for core_atom_id1, core_atom_id2 in it.zip_longest(
        building_block1.get_core_atom_ids(),
        building_block2.get_core_atom_ids(),
    ):
        assert core_atom_id1 == core_atom_id2


def is_clone_building_block(building_block1, building_block2):
    assert building_block1 is not building_block2
    is_equivalent_building_block(building_block1, building_block2)
