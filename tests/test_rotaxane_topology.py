import os
import stk

test_dir = 'rotaxane_topology_tests_output'
if not os.path.exists(test_dir):
    os.mkdir(test_dir)


def test_assembly(amine2, aldehyde2, boronic_acid2, diol2):
    repeat_units = 3

    c1 = stk.Macrocycle([amine2, aldehyde2],
                        stk.Cyclic('AB', [0, 0], repeat_units))

    c2 = stk.Macrocycle([boronic_acid2, diol2],
                        stk.Cyclic('AB', [0, 0], repeat_units))

    axle = stk.Polymer([amine2, aldehyde2],
                       stk.Linear('AB', [0, 0], repeat_units, 'h'))

    r1 = stk.Rotaxane([c1, c2, axle],
                      stk.NRotaxane('ABA', [0, 1, 1], 1))

    path = os.path.join(test_dir, 'rotaxane.mol')
    r1.write(path)

    # No bonds should be made in rotaxane but atoms should be added.

    assert r1.bonds_made == 0

    assert (r1.mol.GetNumAtoms() == axle.mol.GetNumAtoms() +
            c1.mol.GetNumAtoms()*2 + c2.mol.GetNumAtoms())

    assert r1.bb_counter[c1] == 2
    assert r1.bb_counter[c2] == 1
    assert r1.bb_counter[axle] == 1

    assert r1.topology == stk.NRotaxane('ABA', [0, 1, 1], 1)

    assert (r1.mol.GetNumBonds() ==
            c1.mol.GetNumBonds()*2 + c2.mol.GetNumBonds() +
            axle.mol.GetNumBonds())
