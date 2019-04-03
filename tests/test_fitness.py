import stk
import numpy as np


def test_cage(test_mol1):
    assert np.allclose(stk.fitness.cage(test_mol1),
                       [test_mol1.cavity_size(),
                        max(test_mol1.windows()),
                        test_mol1.window_difference(),
                        test_mol1.energy.rdkit()/test_mol1.bonds_made,
                        test_mol1.bb_distortion(),
                        test_mol1.dihedral_strain()],
                       atol=1e-6)


def test_cage_target(test_mol1, amine2):
    fitness = stk.cage_target(
                  test_mol1,
                  amine2,
                  stk.FunctionData('rdkit', forcefield='mmff'),
                  stk.FunctionData('do_not_optimize'))

    expected_fitness = [test_mol1.cavity_size(),
                        test_mol1.window_difference(),
                        test_mol1.bb_distortion(),
                        test_mol1.cavity_size(),
                        test_mol1.window_difference(),
                        test_mol1.bb_distortion(),
                        test_mol1.dihedral_strain()]
    # Binding energy is a point to compare, so don't.
    assert np.allclose(fitness[1:], expected_fitness, atol=1e-6)


def test_cage_c60(test_mol1, c60):
    fitness = stk.cage_c60(
                  test_mol1,
                  c60,
                  stk.FunctionData('rdkit', forcefield='mmff'),
                  stk.FunctionData('do_not_optimize'),
                  n5fold=1,
                  n2fold=1)

    expected_fitness = [
                   -704.24608581,
                   test_mol1.cavity_size(),
                   test_mol1.window_difference(),
                   test_mol1.bb_distortion(),
                   test_mol1.cavity_size(),
                   test_mol1.window_difference(),
                   test_mol1.bb_distortion(),
                   test_mol1.dihedral_strain()]

    assert np.allclose(fitness, expected_fitness, atol=1e-6)
