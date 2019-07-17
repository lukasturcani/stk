import os
import stk
from os.path import join

if not os.path.exists('constructed_molecule_tests_output'):
    os.mkdir('constructed_molecule_tests_output')


def test_caching(amine2, aldehyde2):
    polymer = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2],
        topology=stk.polymer.Linear('AB', [0, 0], 3),
        use_cache=True
    )
    polymer2 = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2],
        topology=stk.polymer.Linear('AB', [0, 0], 3),
        use_cache=True
    )
    assert polymer is polymer2

    polymer3 = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2],
        topology=stk.polymer.Linear('AB', [0, 0], 3),
        use_cache=False
    )
    assert polymer is not polymer3

    polymer4 = stk.Polymer(
        building_blocks=[amine2, aldehyde2],
        topology=stk.polymer.Linear('AB', [1, 0.5], 3),
        use_cache=True
    )
    assert polymer is not polymer4


def test_dump_and_load(tmp_polymer):

    path = join('constructed_molecule_tests_output', 'mol.dump')

    tmp_polymer.test_attr1 = 'something'
    tmp_polymer.test_attr2 = 12
    tmp_polymer.test_attr3 = ['12', 'something', 21]
    tmp_polymer.test_attr4 = 'skip'
    include_attrs = ['test_attr1', 'test_attr2', 'test_attr3']

    # Add some custom atom properties.
    tmp_polymer.atoms[0].some_prop = 'custom atom prop'

    tmp_polymer.dump(path, include_attrs=include_attrs)
    mol2 = stk.Molecule.load(path)

    assert tmp_polymer.__class__ is mol2.__class__
    assert tmp_polymer is not mol2
    assert tmp_polymer.func_groups == mol2.func_groups
    for a1, a2 in zip(tmp_polymer.atoms, mol2.atoms):
        assert a1.__class__ is a2.__class__
        assert vars(a1) == vars(a2)

    for b1, b2 in zip(tmp_polymer.bonds, mol2.bonds):
        assert b1.__class__ is b2.__class__
        assert vars(b1) == vars(b2)

    assert repr(tmp_polymer.topology) == repr(mol2.topology)
    assert tmp_polymer.bonds_made == mol2.bonds_made
    assert tmp_polymer.test_attr1 == mol2.test_attr1
    assert tmp_polymer.test_attr2 == mol2.test_attr2
    assert tmp_polymer.test_attr3 == mol2.test_attr3
    assert not hasattr(mol2, 'test_attr4')
