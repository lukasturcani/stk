from .. import FunctionData, Molecule
from ..ga import cage, cage_target, cage_c60
from os.path import join
import numpy as np

cagemol = Molecule.load(join('data', 'fitness', 'cage.json'))
target = join('data', 'fitness', 'target.pdb')

def test_cage():
    assert np.allclose(cage(cagemol),
        [5.71526484e+00, 3.91281587e+00,
         7.80441714e-04, 1.16913326e+00],
        atol=1e-8)

def test_cage_target():
    assert np.allclose(cage_target(cagemol, target,
        FunctionData('rdkit', forcefield='mmff'),
        FunctionData('rdkit_optimization')),
        [  1.99425829e+02, 7.80441714e-04], atol=1e-8)

def test_cage_c60():
    assert np.allclose(cage_c60(cagemol, target,
            FunctionData('rdkit', forcefield='mmff'),
            FunctionData('rdkit_optimization'), 1, 1),
            [ 1.99443688e+02, 7.80441714e-04], atol=1e-8)
