import os
import stk
import itertools as it
from os.path import join

if not os.path.exists('constructed_molecule_tests_output'):
    os.mkdir('constructed_molecule_tests_output')


def test_init(amine2, aldehyde2):
    polymer = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2],
        topology_graph=stk.polymer.Linear('AB', [0, 0], 3),
        use_cache=True
    )
    polymer2 = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2],
        topology_graph=stk.polymer.Linear('AB', [0, 0], 3),
        use_cache=True
    )
    assert polymer is polymer2

    polymer3 = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2],
        topology_graph=stk.polymer.Linear('AB', [0, 0], 3),
        use_cache=False
    )
    assert polymer is not polymer3

    polymer4 = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2],
        topology_graph=stk.polymer.Linear('AB', [1, 0.5], 3),
        use_cache=True
    )
    assert polymer is not polymer4


def test_dump_and_load(tmp_polymer):

    path = join('constructed_molecule_tests_output', 'mol.dump')

    tmp_polymer.test_attr1 = 'something'
    tmp_polymer.test_attr2 = 12
    tmp_polymer.test_attr3 = ['12', 'something', 21]
    tmp_polymer.test_attr4 = 'skip'

    bb1, bb2 = tmp_polymer.building_block_vertices
    bb1.test_attr1 = 1232
    bb2.test_attr5 = 'alpha'

    include_attrs = [
        'test_attr1',
        'test_attr2',
        'test_attr3',
        'test_attr5'
    ]

    # Add some custom atom properties.
    tmp_polymer.atoms[0].some_prop = 'custom atom prop'
    # Add some custom bond properties
    tmp_polymer.bonds[2].other_prop = 1999

    tmp_polymer.dump(
        path=path,
        include_attrs=include_attrs,
        ignore_missing_attrs=True
    )
    mol2 = stk.Molecule.load(path)

    assert tmp_polymer.__class__ is mol2.__class__
    assert tmp_polymer is not mol2
    fgs = it.zip_longest(tmp_polymer.func_groups, mol2.func_groups)
    for fg1, fg2 in fgs:
        atoms = it.zip_longest(fg1.atoms, fg2.atoms)
        bonders = it.zip_longest(fg1.bonders, fg2.bonders)
        deleters = it.zip_longest(fg1.deleters, fg2.deleters)
        for a1, a2 in it.chain(atoms, bonders, deleters):
            assert a1.__class__ is a2.__class__
            assert a1.id is a1.id

    for a1, a2 in zip(tmp_polymer.atoms, mol2.atoms):
        assert a1.__class__ is a2.__class__
        d1, d2 = vars(a1), vars(a2)
        bb1, bb2 = d1.pop('building_block'), d2.pop('building_block')
        assert d1 == d2
        assert bb1.is_identical(bb2)

    for b1, b2 in zip(tmp_polymer.bonds, mol2.bonds):
        assert b1.__class__ is b2.__class__
        d1, d2 = vars(b1), vars(b2)
        assert repr(d1.pop('atom1')) == repr(d2.pop('atom1'))
        assert repr(d1.pop('atom2')) == repr(d2.pop('atom2'))
        assert d1 == d2

    assert (
        repr(tmp_polymer.topology_graph) == repr(mol2.topology_graph)
    )
    assert tmp_polymer.bonds_made == mol2.bonds_made
    assert tmp_polymer.test_attr1 == mol2.test_attr1
    assert tmp_polymer.test_attr2 == mol2.test_attr2
    assert tmp_polymer.test_attr3 == mol2.test_attr3
    assert not hasattr(mol2, 'test_attr4')

    bbs1 = list(tmp_polymer.building_block_vertices.keys())
    bbs2 = list(mol2.building_block_vertices.keys())
    for bb1, bb2 in it.zip_longest(bbs1, bbs2):
        assert bb1.is_identical(bb2)
        assert bb1 is not bb2
        bb1_count = tmp_polymer.building_block_counter[bb1]
        bb2_count = mol2.building_block_counter[bb2]
        assert bb1_count == bb2_count

    assert bbs2[0].test_attr1 == 1232
    assert bbs2[1].test_attr5 == 'alpha'

    mol3 = stk.Molecule.load(path, use_cache=True)
    assert mol3 is tmp_polymer


def test_is_identical(polymer, tmp_polymer):
    assert polymer is not tmp_polymer
    assert polymer.is_identical(tmp_polymer)
    assert tmp_polymer.is_identical(polymer)
