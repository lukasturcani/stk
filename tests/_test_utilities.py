import itertools as it
import stk
from os.path import join


def _add_test_attrs(mol):
    mol.test_attr1 = 'something'
    mol.test_attr2 = 12
    mol.test_attr3 = ['12', 'something', 21]
    mol.test_attr4 = 'skip'

    bb1, bb2, *_ = mol.get_building_blocks()
    bb1.test_attr1 = 1232
    bb2.test_attr5 = 'alpha'

    # Add some custom atom properties.
    mol.atoms[0].some_prop = 'custom atom prop'
    # Add some custom bond properties
    mol.bonds[2].other_prop = 1999


def _test_func_groups(mol, loaded):
    fgs = it.zip_longest(mol.func_groups, loaded.func_groups)
    for fg1, fg2 in fgs:
        atoms = it.zip_longest(fg1.atoms, fg2.atoms)
        bonders = it.zip_longest(fg1.bonders, fg2.bonders)
        deleters = it.zip_longest(fg1.deleters, fg2.deleters)
        for a1, a2 in it.chain(atoms, bonders, deleters):
            assert a1.__class__ is a2.__class__
            assert a1.id is a1.id


def _test_atoms(mol, loaded):
    for a1, a2 in zip(mol.atoms, loaded.atoms):
        assert a1.__class__ is a2.__class__
        d1, d2 = dict(vars(a1)), dict(vars(a2))
        bb1, bb2 = d1.pop('building_block'), d2.pop('building_block')
        assert d1 == d2
        assert bb1.is_identical(bb2)


def _test_bonds(mol, loaded):
    for b1, b2 in zip(mol.bonds, loaded.bonds):
        assert b1.__class__ is b2.__class__
        d1, d2 = dict(vars(b1)), dict(vars(b2))
        assert repr(d1.pop('atom1')) == repr(d2.pop('atom1'))
        assert repr(d1.pop('atom2')) == repr(d2.pop('atom2'))
        assert d1 == d2


def _test_attrs(mol, loaded):
    assert mol.test_attr1 == loaded.test_attr1
    assert mol.test_attr2 == loaded.test_attr2
    assert mol.test_attr3 == loaded.test_attr3
    assert not hasattr(loaded, 'test_attr4')


def _test_bbs(mol, loaded):
    bbs1 = list(mol.building_block_vertices.keys())
    bbs2 = list(loaded.building_block_vertices.keys())
    for bb1, bb2 in it.zip_longest(bbs1, bbs2):
        assert bb1.is_identical(bb2)
        assert bb1 is not bb2
        bb1_count = mol.building_block_counter[bb1]
        bb2_count = loaded.building_block_counter[bb2]
        assert bb1_count == bb2_count

    assert bbs2[0].test_attr1 == 1232
    assert bbs2[1].test_attr5 == 'alpha'


def _test_dump_and_load(test_dir, mol):
    path = join(
        test_dir, f'{mol.topology_graph.__class__.__name__}.dump'
    )
    _add_test_attrs(mol)
    mol.dump(
        path=path,
        include_attrs=[
            'test_attr1',
            'test_attr2',
            'test_attr3',
            'test_attr5'
        ],
        ignore_missing_attrs=True
    )
    loaded = stk.Molecule.load(path)

    assert mol.__class__ is loaded.__class__
    assert loaded is not mol
    _test_func_groups(mol, loaded)
    _test_atoms(mol, loaded)
    _test_bonds(mol, loaded)

    assert repr(loaded.topology_graph) == repr(mol.topology_graph)
    assert (
        len(loaded.construction_bonds) == len(mol.construction_bonds)
    )
    _test_attrs(mol, loaded)
    _test_bbs(mol, loaded)
    mol3 = stk.Molecule.load(path, use_cache=True)
    mol4 = stk.Molecule.load(path, use_cache=True)
    assert mol3 is mol4
