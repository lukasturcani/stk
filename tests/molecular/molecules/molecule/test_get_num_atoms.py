import stk
import rdkit.Chem.AllChem as rdkit
import pytest


@pytest.mark.parametrize(
    argnames=('molecule', 'num_atoms'),
    argvalues=(
        (stk.BuildingBlock('NCCN'), 12),
    ),
)
def test_get_num_atoms_1(molecule, num_atoms):
    assert molecule.get_num_atoms() == num_atoms


def test_get_num_atoms_2(case_data):
    _test_get_num_atoms_2(case_data.molecule, case_data.smiles)


def _test_get_num_atoms_2(molecule, smiles):
    expected = rdkit.MolFromSmiles(smiles, sanitize=False)
    assert molecule.get_num_atoms() == expected.GetNumAtoms()
