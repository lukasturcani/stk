import os
import stk
import numpy as np
from os.path import join
import itertools as it

test_dir = 'linear_topology_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


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
    assert (
        len(tmp_polymer.construction_bonds) ==
        len(mol2.construction_bonds)
    )
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
    mol4 = stk.Molecule.load(path, use_cache=True)
    assert mol3 is mol4


def test_linear_vertex(tmp_amine2):
    t = stk.polymer.Linear(
        repeating_unit='AB',
        orientations=[0, 0],
        n=3
    )
    v1, v2 = t.vertices[2:4]

    new_coords = v1.place_building_block(tmp_amine2)
    assert np.allclose(
        a=tmp_amine2.get_centroid(tmp_amine2.get_bonder_ids()),
        b=v1.get_position(),
        atol=1e-6
    )
    assert np.allclose(
        a=new_coords,
        b=tmp_amine2.get_position_matrix(),
        atol=1e-6
    )

    new_coords = v2.place_building_block(tmp_amine2)
    assert np.allclose(
        a=tmp_amine2.get_centroid(tmp_amine2.get_bonder_ids()),
        b=v2.get_position(),
        atol=1e-6
    )
    assert np.allclose(
        a=new_coords,
        b=tmp_amine2.get_position_matrix(),
        atol=1e-6
    )


def test_construction(amine2, aldehyde2, boronic_acid2, diol2):
    repeat_units = 3
    monomer_joins = 2*repeat_units - 1

    p1 = stk.ConstructedMolecule(
        building_blocks=[amine2, aldehyde2],
        topology_graph=stk.polymer.Linear(
            repeating_unit='AB',
            orientations=[1, 1],
            n=repeat_units
        )
    )

    p2 = stk.ConstructedMolecule(
        building_blocks=[boronic_acid2, diol2],
        topology_graph=stk.polymer.Linear(
            repeating_unit='AB',
            orientations=[0, 0],
            n=repeat_units
        )
    )

    path = os.path.join(test_dir, 'p1.mol')
    p1.write(path)
    p2.write(path.replace('1', '2'))

    assert len(p1.construction_bonds) == monomer_joins
    assert len(p2.construction_bonds) == monomer_joins*2

    num_monomer_atoms = len(amine2.atoms) + len(aldehyde2.atoms)

    # 3 atoms are lost at each join in p1 due to condensation.
    expected_atoms = num_monomer_atoms*repeat_units - 3*monomer_joins
    assert len(p1.atoms) == expected_atoms

    # 6 atoms are lost at each join due to condensation.
    num_monomer_atoms = len(boronic_acid2.atoms) + len(diol2.atoms)
    expected_atoms = num_monomer_atoms*repeat_units - 6*monomer_joins
    assert len(p2.atoms) == expected_atoms

    assert p1.building_block_counter[amine2] == repeat_units
    assert p1.building_block_counter[aldehyde2] == repeat_units
    assert p2.building_block_counter[boronic_acid2] == repeat_units
    assert p2.building_block_counter[diol2] == repeat_units

    t1 = stk.polymer.Linear('AB', [1, 1], repeat_units)
    assert repr(p1.topology_graph) == repr(t1)
    t2 = stk.polymer.Linear('AB', [0, 0], repeat_units)
    assert repr(p2.topology_graph) == repr(t2)

    expected_bonds = (
        len(amine2.bonds)*repeat_units +
        len(aldehyde2.bonds)*repeat_units -
        monomer_joins*2
    )
    assert len(p1.bonds) == expected_bonds
    expected_bonds = (
        len(boronic_acid2.bonds)*repeat_units +
        len(diol2.bonds)*repeat_units -
        monomer_joins*4
    )
    assert len(p2.bonds) == expected_bonds
