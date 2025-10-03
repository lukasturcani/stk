import pytest

import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        lambda name: CaseData(
            molecules=[
                stk.BuildingBlock(smiles="CC"),
                stk.BuildingBlock(smiles="CNC"),
                stk.BuildingBlock(smiles="CNNC"),
                stk.BuildingBlock(smiles="CNC1CCC1C"),
            ],
            property_dicts=[
                {"1": 2, "3": "astr"},
                {"1": 2, "3": "bstr"},
                {"1": 2, "3": "cstr"},
                {"1": 3, "3": "cstr", "2": {"h": "w"}},
            ],
            expected_count=4,
            name=name,
        ),
    )
)
def molecule(request: pytest.FixtureRequest) -> CaseData:
    return request.param(
        f"{request.fixturename}{request.param_index}",
    )
