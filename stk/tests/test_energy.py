from os.path import join
from ..ga import Population
from ..convenience_tools import FunctionData
from ..molecular import Molecule
import numpy as np


pop = Population.load(join('data', 'macromolecule', 'mm.json'),
                      Molecule.from_dict)
mol, mol2 = pop[:2]


def test_rdkit():
    assert np.isclose(mol.energy.rdkit('uff'), -0.84, atol=0.01)
    assert np.isclose(mol.energy.rdkit('mmff'), 3476.86, atol=0.02)


def test_formation():
    assert np.isclose(
            mol.energy.formation(
                        FunctionData('rdkit', forcefield='uff'),
                        [(2, mol2)]),
            0,
            atol=18.318412936376188)


def test_pseudoformation():
    assert np.isclose(
            mol.energy.pseudoformation(
                    FunctionData('rdkit', forcefield='uff')),
            0,
            atol=17.75271336716965)


def test_logging():
    assert len(mol.energy.values) != 0
