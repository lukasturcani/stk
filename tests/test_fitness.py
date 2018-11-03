import stk
from os.path import join
import numpy as np
import pytest

cagemol = stk.Molecule.load(join('data', 'fitness', 'cage.json'))
target = join('data', 'fitness', 'target.pdb')


def test_cage():
    assert isinstance(stk.fitness.cage(cagemol), np.ndarray)


def test_cage_target():
    # Because no optimization function is used this function will fail.
    with pytest.raises(ValueError) as ex:
        isinstance(stk.cage_target(cagemol, target,
                   stk.FunctionData('rdkit', forcefield='mmff'),
                   stk.FunctionData('do_not_optimize')), np.ndarray)
        # The 3rd fitness parameter should be the only one that failed.
        ex.args[1].pop(2)
        assert all(x is not None for x in ex.args[1])


def test_cage_c60():
    # Because no optimization function is used this function will fail.
    with pytest.raises(ValueError) as ex:
        isinstance(stk.cage_c60(cagemol, target,
                   stk.FunctionData('rdkit', forcefield='mmff'),
                   stk.FunctionData('do_not_optimize'), 1, 1),
                   np.ndarray)
        # The 3rd fitness parameter should be the only one that failed.
        ex.args[1].pop(2)
        assert all(x is not None for x in ex.args[1])
