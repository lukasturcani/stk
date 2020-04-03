import pytest
import numpy as np
import stk

from ..case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            fitness_calculator=stk.PropertyVector(
                stk.Molecule.get_num_atoms,
                stk.Molecule.get_num_bonds,
                stk.Molecule.get_maximum_diameter,
            ),
            molecule=stk.BuildingBlock('BrCCBr').with_position_matrix(
                position_matrix=np.array([
                    [0, 0, 0],
                    [10, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                ], dtype=np.float64),
            ),
            fitness_value=(8, 7, 10),
        ),
    ),
)
def property_vector(request):
    return request.param
