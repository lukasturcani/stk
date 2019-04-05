import os
import stk

if not os.path.exists('macromolecule_tests_output'):
    os.mkdir('macromolecule_tests_output')


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


def test_comparison():
    """
    Checks ``==``, ``>``, ``>=``, etc. operators.

    """

    a = stk.MacroMolecule.__new__(stk.MacroMolecule)
    a.fitness = 1

    b = stk.MacroMolecule.__new__(stk.MacroMolecule)
    b.fitness = 1

    c = stk.MacroMolecule.__new__(stk.MacroMolecule)
    c.fitness = 2

    # Comparison operators should compare their fitness.
    assert not a < b
    assert a <= b
    assert a == b
    assert c > b
    assert c >= a


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


def test_json_init(polymer):
    path = os.path.join('macromolecule_tests_output', 'mol.json')
    polymer.dump(path)
    mol2 = stk.Molecule.load(path, stk.Molecule.from_dict)

    assert polymer is not mol2
    assert polymer.atom_props == mol2.atom_props
    assert polymer.bb_counter == mol2.bb_counter
    assert polymer.topology == mol2.topology
    assert polymer.bonds_made == mol2.bonds_made
    assert polymer.unscaled_fitness == mol2.unscaled_fitness
    assert polymer.progress_params == mol2.progress_params
    assert all(bb1.same(bb2) and
               bb1.func_groups == bb2.func_groups
               for bb1, bb2 in
               zip(polymer.building_blocks, mol2.building_blocks))
