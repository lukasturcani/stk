import stk
import pytest
import rdkit.Chem.AllChem as rdkit


@pytest.mark.parametrize(
    argnames=('molecule', 'num_bonds'),
    argvalues=(
        (stk.BuildingBlock('NCCN'), 11),
    ),
)
def test_get_num_bonds(molecule, num_bonds):
    assert molecule.get_num_bonds() == num_bonds


def test_get_num_bonds_2(case_data):
    _test_get_num_bonds_2(case_data.molecule, case_data.smiles)


def _test_get_num_bonds_2(molecule, smiles):
    expected = rdkit.AddHs(rdkit.MolFromSmiles(smiles))
    assert molecule.get_num_bonds() == expected.GetNumBonds()
