import os
import stk

test_dir = 'linear_topology_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def test_constrution(amine2, aldehyde2, boronic_acid2, diol2):
    repeat_units = 3
    monomer_joins = 2*repeat_units - 1

    p1 = stk.Polymer([amine2, aldehyde2],
                     stk.Linear('AB', [0, 0], repeat_units, 'h'))

    p2 = stk.Polymer([amine2, aldehyde2],
                     stk.Linear('AB', [1, 1], repeat_units, 'fg'))

    p3 = stk.Polymer([boronic_acid2, diol2],
                     stk.Linear('AB', [0, 0], repeat_units, 'h'))

    path = os.path.join(test_dir, 'p1.mol')
    p1.write(path)
    p2.write(path.replace('1', '2'))
    p3.write(path.replace('1', '3'))

    assert p1.bonds_made == monomer_joins
    assert p2.bonds_made == monomer_joins
    assert p3.bonds_made == monomer_joins*2

    monomer_atom_count = (amine2.mol.GetNumAtoms() +
                          aldehyde2.mol.GetNumAtoms())

    # 3 atoms are lost at each join in p1 and p2 due to condensation.
    # In p1 an atom is gained due to replacing aldehyde fg with Hs.
    assert (p1.mol.GetNumAtoms() ==
            monomer_atom_count*repeat_units - 3*monomer_joins + 1)

    assert (p2.mol.GetNumAtoms() ==
            monomer_atom_count*repeat_units - 3*monomer_joins)

    # p3 loses two atoms because boronic acid OH groups are repalced
    # with H.
    monomer_atom_count = (boronic_acid2.mol.GetNumAtoms() +
                          diol2.mol.GetNumAtoms())

    assert (p3.mol.GetNumAtoms() ==
            monomer_atom_count*repeat_units - 6*monomer_joins - 2)

    assert p1.bb_counter[amine2] == repeat_units
    assert p1.bb_counter[aldehyde2] == repeat_units
    assert p2.bb_counter[amine2] == repeat_units
    assert p2.bb_counter[aldehyde2] == repeat_units
    assert p3.bb_counter[boronic_acid2] == repeat_units
    assert p3.bb_counter[diol2] == repeat_units

    assert p1.topology == stk.Linear('AB', [0, 0], repeat_units, 'h')
    assert p2.topology == stk.Linear('AB', [1, 1], repeat_units, 'fg')
    assert p3.topology == stk.Linear('AB', [0, 0], repeat_units, 'h')

    assert (p1.mol.GetNumBonds() ==
            amine2.mol.GetNumBonds()*repeat_units +
            aldehyde2.mol.GetNumBonds()*repeat_units -
            monomer_joins*2 + 1)
    assert (p2.mol.GetNumBonds() ==
            amine2.mol.GetNumBonds()*repeat_units +
            aldehyde2.mol.GetNumBonds()*repeat_units -
            monomer_joins*2)
    assert (p3.mol.GetNumBonds() ==
            boronic_acid2.mol.GetNumBonds()*repeat_units +
            diol2.mol.GetNumBonds()*repeat_units -
            monomer_joins*4 - 2)
