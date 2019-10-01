import os
import stk
from collections import Counter
import itertools as it


if not os.path.exists('constructed_molecule_tests_output'):
    os.mkdir('constructed_molecule_tests_output')


def test_init(amine2, aldehyde2):
    polymer = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2],
        topology_graph=stk.polymer.Linear('AB', 3),
        use_cache=True
    )
    polymer2 = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2],
        topology_graph=stk.polymer.Linear('AB', 3),
        use_cache=True
    )
    assert polymer is polymer2

    polymer3 = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2],
        topology_graph=stk.polymer.Linear('AB', 3),
        use_cache=False
    )
    assert polymer is not polymer3

    polymer4 = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2],
        topology_graph=stk.polymer.Linear('AB', 3, (1, 0.5)),
        use_cache=True
    )
    assert polymer is not polymer4


def test_get_identity_key(
    polymer,
    tmp_polymer,
    amine2,
    aldehyde3,
    amine2_alt1
):
    assert polymer is not tmp_polymer
    assert polymer.get_identity_key() == tmp_polymer.get_identity_key()

    four_plus_six = stk.cage.FourPlusSix()
    cage1 = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde3, amine2_alt1],
        topology_graph=four_plus_six,
        building_block_vertices={
            aldehyde3: four_plus_six.vertices[:4],
            amine2: four_plus_six.vertices[4:6],
            amine2_alt1: four_plus_six.vertices[6:]
        }
    )
    cage2 = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde3, amine2_alt1],
        topology_graph=four_plus_six,
        building_block_vertices={
            aldehyde3: four_plus_six.vertices[:4],
            amine2: four_plus_six.vertices[4:5],
            amine2_alt1: four_plus_six.vertices[5:]
        }
    )
    assert cage1.get_identity_key() != cage2.get_identity_key()


def _test_atoms(atoms):
    for a1, a2 in atoms:
        assert a1 is not a2
        assert a1.id == a2.id
        assert a1.__class__ is a2.__class__
        assert a1.charge == a2.charge


def _test_bonds(bonds):
    for b1, b2 in bonds:
        assert b1 is not b2
        assert b1.__class__ is b2.__class__
        assert b1.atom1.id == b2.atom1.id
        assert b1.order == b2.order
        assert b1.periodicity == b2.periodicity


def test_clone(polymer):
    clone = polymer.clone()

    atoms = it.zip_longest(polymer.atoms, clone.atoms)
    _test_atoms(atoms)

    bonds = it.zip_longest(polymer.bonds, clone.bonds)
    _test_bonds(bonds)

    fgs = it.zip_longest(polymer.func_groups, clone.func_groups)
    for fg1, fg2 in fgs:
        atoms = it.zip_longest(fg1.atoms, fg2.atoms)
        _test_atoms(atoms)
        bonders = it.zip_longest(fg1.bonders, fg2.bonders)
        _test_atoms(bonders)
        deleters = it.zip_longest(fg1.deleters, fg2.deleters)
        _test_atoms(deleters)

    assert repr(polymer.topology_graph) == repr(clone.topology_graph)
    bb_vertices = it.zip_longest(
        polymer.building_block_vertices.items(),
        clone.building_block_vertices.items(),
    )
    for (bb1, vertices1), (bb2, vertices2) in bb_vertices:
        assert bb1.get_identity_key() == bb2.get_identity_key()
        for v1, v2 in it.zip_longest(vertices1, vertices2):
            assert v1.id == v2.id

    assert (
        polymer.building_block_counter == clone.building_block_counter
    )
    construction_bonds = it.zip_longest(
        polymer.construction_bonds,
        clone.construction_bonds,
    )
    _test_bonds(construction_bonds)


def test_get_building_blocks(
    amine2,
    amine2_alt1,
    aldehyde2,
    aldehyde2_alt1,
    aldehyde3,
    polymer,
    polymer_alt1,
    four_plus_six
):
    polymer_bbs = Counter(polymer.get_building_blocks())
    assert len(polymer_bbs) == 2
    assert polymer_bbs[amine2] == 1
    assert polymer_bbs[aldehyde2] == 1

    polymer_alt1_bbs = Counter(polymer_alt1.get_building_blocks())
    assert len(polymer_alt1_bbs) == 2
    assert polymer_alt1_bbs[amine2_alt1] == 1
    assert polymer_alt1_bbs[aldehyde2_alt1] == 1

    four_plus_six_bbs = Counter(four_plus_six.get_building_blocks())
    assert len(four_plus_six_bbs) == 2
    assert four_plus_six_bbs[amine2] == 1
    assert four_plus_six_bbs[aldehyde3] == 1
