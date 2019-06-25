import os
import stk
import pickle
from os.path import join

if not os.path.exists('constructed_molecule_tests_output'):
    os.mkdir('constructed_molecule_tests_output')


def test_building_block_cores(polymer):
    # Check that the yielded rdkit molecules match the cores of the
    # building block molecules.
    for i in range(len(polymer.building_blocks)):
        for frag in polymer.building_block_cores(i):
            bb1match = len(frag.GetSubstructMatch(
                           polymer.building_blocks[0].core()))
            bb2match = len(frag.GetSubstructMatch(
                           polymer.building_blocks[1].core()))
            nfrag_atoms = frag.GetNumAtoms()
            assert bb1match == nfrag_atoms or bb2match == nfrag_atoms


def test_bb_distortion(polymer):
    assert isinstance(polymer.bb_distortion(), float)


def test_caching(amine2, aldehyde2):
    # Other tests assume that the cache is turned off. Make sure
    # that this test restores cache to off after it finishes.
    try:
        stk.OPTIONS['cache'] = True
        polymer = stk.Polymer([amine2, aldehyde2],
                              stk.Linear('AB', [0, 0], 3))
        polymer2 = stk.Polymer([amine2, aldehyde2],
                               stk.Linear('AB', [0, 0], 3))
        assert polymer is polymer2

        stk.OPTIONS['cache'] = False
        polymer3 = stk.Polymer([amine2, aldehyde2],
                               stk.Linear('AB', [1, 0.5], 3))
        assert polymer is not polymer3
    except Exception:
        raise
    finally:
        stk.OPTIONS['cache'] = False


def test_save_rdkit_atom_props(tmp_amine2):
    assert 'int_test' not in tmp_amine2.atom_props[0]
    assert 'int_test2' not in tmp_amine2.atom_props[0]
    assert 'str_test' not in tmp_amine2.atom_props[2]

    tmp_amine2.mol.GetAtomWithIdx(0).SetIntProp('int_test', 121)
    tmp_amine2.mol.GetAtomWithIdx(2).SetProp('str_test', 'value')
    tmp_amine2.save_rdkit_atom_props({'int_test'})

    assert tmp_amine2.atom_props[0]['int_test'] == 121
    assert 'int_test2' not in tmp_amine2.atom_props[0]
    assert 'str_test' not in tmp_amine2.atom_props[2]

    tmp_amine2.save_rdkit_atom_props({'int_test2', 'str_test'})
    assert tmp_amine2.atom_props[0]['int_test'] == 121
    assert 'int_test2' not in tmp_amine2.atom_props[0]
    assert tmp_amine2.atom_props[2]['str_test'] == 'value'


def test_json_init(tmp_polymer):

    path = join('constructed_molecule_tests_output', 'mol.json')

    tmp_polymer.test_attr1 = 'something'
    tmp_polymer.test_attr2 = 12
    tmp_polymer.test_attr3 = ['12', 'something', 21]
    include_attrs = ['test_attr1', 'test_attr2', 'test_attr3']

    tmp_polymer.dump(path, include_attrs=include_attrs)
    mol2 = stk.Molecule.load(path, stk.Molecule.from_dict)

    assert tmp_polymer.__class__ is mol2.__class__
    assert tmp_polymer is not mol2
    assert tmp_polymer.func_groups == mol2.func_groups
    assert tmp_polymer.atom_props == mol2.atom_props

    matches = 0
    for key1, value1 in tmp_polymer.bb_counter.items():
        for key2, value2 in mol2.bb_counter.items():
            if key1.same(key2):
                assert value1 == value2
                matches += 1
    assert len(tmp_polymer.bb_counter) == len(mol2.bb_counter)
    assert len(mol2.bb_counter) == matches

    assert tmp_polymer.topology == mol2.topology
    assert tmp_polymer.bonds_made == mol2.bonds_made
    assert all(bb1.same(bb2) and
               bb1.func_groups == bb2.func_groups
               for bb1, bb2 in
               zip(tmp_polymer.building_blocks, mol2.building_blocks))
    assert tmp_polymer.test_attr1 == mol2.test_attr1
    assert tmp_polymer.test_attr2 == mol2.test_attr2
    assert tmp_polymer.test_attr3 == mol2.test_attr3


def test_pickle(polymer):
    result = pickle.loads(pickle.dumps(polymer))
    assert result.same(polymer)
