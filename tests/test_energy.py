from os.path import join
from ..molecular import Molecule
import numpy as np

path = join('data', 'cage', 'cage.json')
with open(path, 'r') as f:
    mol = Molecule.load(f.read())

def test_rdkit():
    assert np.isclose(mol.energy.rdkit('uff'), -0.84, atol=0.01)
    assert np.isclose(mol.energy.rdkit('mmff'), 3476.86, atol=0.02)

def test_formation():
    ...

def test_pseudoformation():
    ...

def test_logging():
    assert len(mol.energy.values.keys()) == 2
