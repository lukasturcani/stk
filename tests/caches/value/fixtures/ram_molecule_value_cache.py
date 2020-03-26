import pytest
import stk

from ..case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            cache=stk.RamMoleculeValueCache(),
            molecule=stk.BuildingBlock('BrCCBr'),
            value=12,
        ),
    ),
)
def ram_molecule_value_cache(request):
    return request.param
