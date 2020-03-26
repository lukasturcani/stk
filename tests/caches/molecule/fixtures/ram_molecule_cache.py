import pytest
import rdkit.Chem.AllChem as rdkit
import stk

from ..case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            cache=stk.RamMoleculeCache(),
            molecule=stk.BuildingBlock('BrCCBr'),
            key=rdkit.MolToInchiKey(rdkit.MolFromSmiles('BrCCBr')),
        ),
    ),
)
def ram_molecule_cache(request):
    return request.param
