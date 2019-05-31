import os
import stk

test_dir = 'cyclic_topology_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def test_assembly(amine2, aldehyde2, boronic_acid2, diol2):
    repeat_units = 3
    monomer_joins = 2*repeat_units

    c1 = stk.Macrocycle([amine2, aldehyde2],
                        stk.Cyclic('AB', [0, 0], repeat_units))

    c2 = stk.Macrocycle([boronic_acid2, diol2],
                        stk.Cyclic('AB', [0, 0], repeat_units))

    path = os.path.join(test_dir, 'c1.mol')
    c1.write(path)
    c2.write(path.replace('1', '2'))

    assert c1.bonds_made == monomer_joins
    assert c2.bonds_made == monomer_joins*2

    monomer_atom_count = (amine2.mol.GetNumAtoms() +
                          aldehyde2.mol.GetNumAtoms())

    # 3 atoms are lost at each join in c1 and c2 due to condensation.
    assert (c1.mol.GetNumAtoms() ==
            monomer_atom_count*repeat_units - 3*monomer_joins)

    # c2 loses 2 atoms because OH groups are replaced with H.
    monomer_atom_count = (boronic_acid2.mol.GetNumAtoms() +
                          diol2.mol.GetNumAtoms())

    assert (c2.mol.GetNumAtoms() ==
            monomer_atom_count*repeat_units - 6*monomer_joins)

    assert c1.bb_counter[amine2] == repeat_units
    assert c1.bb_counter[aldehyde2] == repeat_units
    assert c2.bb_counter[boronic_acid2] == repeat_units
    assert c2.bb_counter[diol2] == repeat_units

    assert c1.topology == stk.Cyclic('AB', [0, 0], repeat_units)
    assert c2.topology == stk.Cyclic('AB', [0, 0], repeat_units)

    assert (c1.mol.GetNumBonds() ==
            amine2.mol.GetNumBonds()*repeat_units +
            aldehyde2.mol.GetNumBonds()*repeat_units -
            monomer_joins*2)
    assert (c2.mol.GetNumBonds() ==
            boronic_acid2.mol.GetNumBonds()*repeat_units +
            diol2.mol.GetNumBonds()*repeat_units -
            monomer_joins*4)
