from os.path import join
from ..molecular import Molecule
import numpy as np
import json

path = join('data', 'macromolecule', 'mm.json')
with open(path, 'r') as f:
    mol = Molecule.load(json.load(f))

def test_rdkit():
    print(mol)
    assert np.isclose(mol.energy.rdkit('uff'), -0.84, atol=0.01)
    assert np.isclose(mol.energy.rdkit('mmff'), 3476.86, atol=0.02)

def test_formation():
    ...

def test_pseudoformation():
    ...

def test_logging():
    print(mol)
    assert len(mol.energy.values.keys()) == 2
