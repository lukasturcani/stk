from os.path import join
from ..ga import Population
import numpy as np


mol = Population.load(join('data', 'macromolecule', 'mm.json'))[0]

def test_rdkit():
    print(mol)
    assert np.isclose(mol.energy.rdkit('uff'), -0.84, atol=0.01)
    assert np.isclose(mol.energy.rdkit('mmff'), 3476.86, atol=0.02)

def test_formation():
    ...

def test_pseudoformation():
    ...

def test_logging():
    assert len(mol.energy.values.keys()) == 2
