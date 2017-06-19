from .. import FunctionData, Molecule
from ..ga.fitness import cage, cage_target, cage_c60
from os.path import join
import numpy as np
import pytest

cagemol = Molecule.load(join('data', 'fitness', 'cage.json'))
target = join('data', 'fitness', 'target.pdb')


def test_cage():
    assert isinstance(cage(cagemol), np.ndarray)


def test_cage_target():
    # Because no optimization function is used this function will fail.
    with pytest.raises(ValueError) as ex:
        isinstance(cage_target(cagemol, target,
                   FunctionData('rdkit', forcefield='mmff'),
                   FunctionData('do_not_optimize')), np.ndarray)
        # The 3rd fitness parameter should be the only one that failed.
        ex.args[1].pop(2)
        assert all(x is not None for x in ex.args[1])


def test_cage_c60():
    # Because no optimization function is used this function will fail.
    with pytest.raises(ValueError) as ex:
        isinstance(cage_c60(cagemol, target,
                   FunctionData('rdkit', forcefield='mmff'),
                   FunctionData('do_not_optimize'), 1, 1),
                   np.ndarray)
        # The 3rd fitness parameter should be the only one that failed.
        ex.args[1].pop(2)
        assert all(x is not None for x in ex.args[1])
