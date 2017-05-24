from .. import FunctionData, Molecule
from ..ga.fitness import cage, cage_target, cage_c60
from os.path import join
import numpy as np

cagemol = Molecule.load(join('data', 'fitness', 'cage.json'))
target = join('data', 'fitness', 'target.pdb')


def test_cage():
    assert isinstance(cage(cagemol), np.ndarray)


def test_cage_target():
    assert isinstance(cage_target(cagemol, target,
                      FunctionData('rdkit', forcefield='mmff'),
                      FunctionData('do_not_optimize')), np.ndarray)


def test_cage_c60():
    assert isinstance(cage_c60(cagemol, target,
                      FunctionData('rdkit', forcefield='mmff'),
                      FunctionData('do_not_optimize'), 1, 1),
                      np.ndarray)
