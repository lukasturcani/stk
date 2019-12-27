import itertools as it

from ..utilities import is_equivalent_atom, is_equivalent_molecule


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

    bonders = it.zip_longest(
        functional_group1.get_bonder_ids(),
        functional_group2.get_bonder_ids(),
    )
    for bonder1, bonder2 in bonders:
        is_equivalent_atom(bonder1, bonder2)

    deleters = it.zip_longest(
        functional_group1.get_deleter_ids(),
        functional_group2.get_deleter_ids(),
    )
    for deleter1, deleter2 in deleters:
        is_equivalent_atom(deleter1, deleter2)


def is_equivalent_building_block(building_block1, building_block2):
    is_equivalent_molecule(building_block1, building_block2)
    functional_groups = it.zip_longest(
        building_block1.get_functional_groups(),
        building_block2.get_functional_groups(),
    )
    for functional_group1, functional_group2 in functional_groups:
        is_equivalent_functional_group(
            functional_group1=functional_group1,
            functional_group2=functional_group2,
        )
