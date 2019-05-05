import stk
import numpy as np
import pytest


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


def test_cache_use(tmp_amine2):
    calc = stk.PropertyVector(lambda mol, conformer: 1,
                              use_cache=False)
    calc.fitness(tmp_amine2)

    # Since use_cache is False the cache should be empty.
    assert not calc.cache

    # To test that the cache is not being used, put a random object
    # into it, and test that it was not returned.
    obj = object()
    calc.cache[(tmp_amine2.key, 1)] = obj
    assert calc.fitness(tmp_amine2, 1) is not obj

    # Test that the cache is being filled when use_cache is True.
    calc = stk.PropertyVector(lambda mol, conformer: 1, use_cache=True)
    assert not calc.cache
    calc.fitness(tmp_amine2)
    assert calc.cache

    # Test that the cache is being used by putting a random object into
    # it and making sure it gets returned.
    calc.cache[(tmp_amine2.key, 1)] = obj
    assert calc.fitness(tmp_amine2, 1) is obj


def test_attribute_creation(tmp_amine2):
    calc = stk.PropertyVector(lambda mol, conformer: 1,
                              use_cache=False)
    assert not hasattr(tmp_amine2, 'fitness')
    calc.fitness(tmp_amine2)
    assert tmp_amine2.fitness == [1]

    calc = stk.PropertyVector(lambda mol, conformer: 2)
    del tmp_amine2.fitness
    assert not hasattr(tmp_amine2, 'fitness')
    calc.fitness(tmp_amine2)
    assert tmp_amine2.fitness == [2]


def test_raising_fitness_calculator(fitness_calculator, tmp_amine2):
    never_raiser = stk.RaisingFitnessCalculator(fitness_calculator, 0)
    never_raiser.fitness(tmp_amine2)

    always_raiser = stk.RaisingFitnessCalculator(fitness_calculator, 1)
    with pytest.raises(stk.RaisingFitnessCalculatorError):
        always_raiser.fitness(tmp_amine2)
