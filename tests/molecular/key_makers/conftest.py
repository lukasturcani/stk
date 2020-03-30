import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            key_maker=stk.Inchi(),
            molecule=stk.BuildingBlock('NCCN'),
            key_name='InChI',
            key='InChI=1S/C2H8N2/c3-1-2-4/h1-4H2',
        ),
        CaseData(
            key_maker=stk.InchiKey(),
            molecule=stk.BuildingBlock('NCCN'),
            key_name='InChIKey',
            key='PIICEJLVQHRZGT-UHFFFAOYSA-N',
        ),
        CaseData(
            key_maker=stk.MoleculeKeyMaker(
                key_name='NumAtoms',
                get_key=lambda molecule: molecule.get_num_atoms(),
            ),
            molecule=stk.BuildingBlock('NCCN'),
            key_name='NumAtoms',
            key=12,
        ),
    ),
)
def case_data(request):
    return request.param
