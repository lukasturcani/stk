import os
import stk

if not os.path.exists('macromolecule_tests_output'):
    os.mkdir('macromolecule_tests_output')

bb1 = stk.StructUnit2.smiles_init('Nc1ccc(N)cc1', 'amine')
bb2 = stk.StructUnit2.smiles_init('O=Cc1cc2ccc3cc(C=O)cc4ccc(c1)c2c34',
                                  'aldehyde')
mol = stk.Polymer([bb1, bb2], stk.Linear('AB', [0.5, 0.5], 3))


def test_building_block_cores():
    # Check that the yielded rdkit molecules match the cores of the
    # building block molecules.
    for i in range(len(mol.building_blocks)):
        for frag in mol.building_block_cores(i):
            bb1match = len(frag.GetSubstructMatch(
                           mol.building_blocks[0].core()))
            bb2match = len(frag.GetSubstructMatch(
                           mol.building_blocks[1].core()))
            nfrag_atoms = frag.GetNumAtoms()
            assert bb1match == nfrag_atoms or bb2match == nfrag_atoms


def test_bb_distortion():
    assert isinstance(mol.bb_distortion(), float)


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


def test_caching():
    mol2 = stk.Polymer([bb2, bb1], stk.Linear('AB', [0.5, 0.5], 3))
    assert mol is mol2

    mol3 = stk.Polymer([bb1, bb2], stk.Linear('AB', [1, 0.5], 3))
    assert mol is not mol3


def test_json_init():
    try:
        path = os.path.join('macromolecule_tests_output', 'mol.json')
        mol.dump(path)
        stk.CACHE_SETTINGS['ON'] = False
        mol2 = stk.Molecule.load(path, stk.Molecule.from_dict)
        stk.CACHE_SETTINGS['ON'] = True

        assert mol is not mol2
        assert mol.bonder_ids == mol2.bonder_ids
        assert mol.atom_props == mol2.atom_props
        assert mol.bb_counter == mol2.bb_counter
        assert mol.topology == mol2.topology
        assert mol.bonds_made == mol2.bonds_made
        assert mol.unscaled_fitness == mol2.unscaled_fitness
        assert mol.progress_params == mol2.progress_params
        assert all(bb1.same(bb2) and bb1.func_grp.name == bb2.func_grp.name
                   for bb1, bb2 in
                   zip(mol.building_blocks, mol2.building_blocks))
    except Exception:
        raise
    finally:
        stk.CACHE_SETTINGS['ON'] = True
