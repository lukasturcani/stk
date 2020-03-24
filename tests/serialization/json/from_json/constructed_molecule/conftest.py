import pytest
import stk
import numpy as np

from .case_data import CaseData
"""

@pytest.fixture(
    params=(
        CaseData(
            molecule_json={
            },
            constructed_molecule_json={
            },
            position_matrix=np.array([
            ]),
            building_blocks=(
            ),
            molecule=stk.ConstructedMolecule(
            ),
        ),
    ),
)
def case_data(request):
    return request.param
"""
