import stk
import numpy as np


def test_property_vector(test_mol1):
    def window_variance(mol, conformer):
        return mol.window_variance(conformer)

    def cavity_size(mol, conformer):
        return mol.cavity_size(conformer)

    def bb_distortion(mol, conformer):
        return mol.bb_distortion(conformer)

    def dihedral_strain(mol, conformer):
        return mol.dihedral_strain(conformer)

    fitness_calculator = stk.PropertyVector(window_variance,
                                            cavity_size,
                                            bb_distortion,
                                            dihedral_strain)
    assert np.allclose(fitness_calculator.fitness(test_mol1),
                       [test_mol1.window_variance(),
                        test_mol1.cavity_size(),
                        test_mol1.bb_distortion(),
                        test_mol1.dihedral_strain()],
                       atol=1e-6)
