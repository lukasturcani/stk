import os
import stk

if not os.path.exists('linear_topology_tests'):
    os.mkdir('linear_topology_tests')


def test_assembly(mol5, mol6, mol7, mol8):
    p1 = stk.Polymer([mol5, mol6], stk.Linear('AB', [0, 0], 3))
    p2 = stk.Polymer([mol5, mol6], stk.Linear('AB', [1, 1], 3, 'fg'))
    p3 = stk.Polymer([mol7, mol8], stk.Linear('AB', [0, 0], 3))

    path = os.path.join('linear_topology_tests', 'p1.mol')
    p1.write(path)
    p2.write(path.replace('1', '2'))
    p3.write(path.replace('1', '3'))
