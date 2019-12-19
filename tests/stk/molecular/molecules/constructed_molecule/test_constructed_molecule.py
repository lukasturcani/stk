import stk
import itertools as it
import pytest


def is_equivalent_atom(atom1, atom2):
    return (
        atom1.id == atom2.id
        and atom1.charge == atom2.charge
        and atom1.__class__ is atom2.__class__
    )


def is_equivalent_bond(bond1, bond2):
    return (
        bond1.__class__ is bond2.__class__
        and bond1.order == bond2.order
        and is_equivalent_atom(bond1.atom1, bond2.atom1)
        and is_equivalent_atom(bond1.atom2, bond2.atom2)
        and bond1.periodicity == bond2.periodicity
    )


def is_equivalent_fg(fg1, fg2):
    equivalent_atoms = all(
        id1 == id2
        for id1, id2
        in it.zip_longest(fg1.get_atom_ids(), fg2.get_atom_ids())
    )
    equivalent_bonders = all(
        id1 == id2
        for id1, id2
        in it.zip_longest(fg1.get_bonder_ids(), fg2.get_bonder_ids())
    )
    equivalent_deleters = all(
        id1 == id2
        for id1, id2
        in it.zip_longest(fg1.get_deleter_ids(), fg2.get_deleter_ids())
    )
    return (
        fg1.fg_type is fg2.fg_type
        and equivalent_atoms
        and equivalent_bonders
        and equivalent_deleters
    )


def test_clone(constructed_molecule):
    clone = constructed_molecule.clone()
    assert (
        repr(constructed_molecule.topology_graph)
        == repr(clone.topology_graph)
    )

    fgs = it.zip_longest(
        clone.func_groups,
        constructed_molecule.func_groups,
    )
    for fg1, fg2 in fgs:
        assert is_equivalent_fg(fg1, fg2)

    assert (
        clone.building_block_counter
        == constructed_molecule.building_block_counter
    )

    construction_bonds = it.zip_longest(
        clone.construction_bonds,
        constructed_molecule.construction_bonds,
    )
    for bond1, bond2 in construction_bonds:
        assert is_equivalent_bond(bond1, bond2)


class TestGetBuildingBlock:
    def case1():
        building_blocks = [stk.BuildingBlock('BrCCBr', ['bromine'])]
        constructed_molecule = stk.ConstructedMolecule(
            building_blocks=building_blocks,
            topology_graph=stk.polymer.Linear('A', 3),
        )
        return constructed_molecule, building_blocks

    @pytest.mark.parametrize(
        argnames=(
            'constructed_molecule',
            'expected_building_blocks',
        ),
        argvalues=(
            case1(),
        ),
    )
    def test(
        self,
        constructed_molecule,
        expected_building_blocks,
    ):
        building_blocks = it.zip_longest(
            constructed_molecule.get_building_blocks(),
            expected_building_blocks,
        )
        for bb1, bb2 in building_blocks:
            assert bb1 is bb2
