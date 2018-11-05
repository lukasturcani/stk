import stk
import numpy as np
import pytest


def test_cage(cage):
    assert isinstance(stk.fitness.cage(cage), np.ndarray)


def test_cage_target(cage, target):
    # Because no optimization function is used this function will fail.
    with pytest.raises(ValueError) as ex:
        isinstance(stk.cage_target(cage, target,
                   stk.FunctionData('rdkit', forcefield='mmff'),
                   stk.FunctionData('do_not_optimize')), np.ndarray)
        # The 3rd fitness parameter should be the only one that failed.
        ex.args[1].pop(2)
        assert all(x is not None for x in ex.args[1])


def test_cage_c60(cage, target):
    # Because no optimization function is used this function will fail.
    with pytest.raises(ValueError) as ex:
        isinstance(stk.cage_c60(cage, target,
                   stk.FunctionData('rdkit', forcefield='mmff'),
                   stk.FunctionData('do_not_optimize'), 1, 1),
                   np.ndarray)
        # The 3rd fitness parameter should be the only one that failed.
        ex.args[1].pop(2)
        assert all(x is not None for x in ex.args[1])
