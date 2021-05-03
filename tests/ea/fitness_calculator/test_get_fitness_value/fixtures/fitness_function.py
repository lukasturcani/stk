import pytest
import stk

from ..case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            fitness_calculator=stk.FitnessFunction(
                fitness_function=stk.Molecule.get_num_atoms,
            ),
            molecule=stk.BuildingBlock('BrCCBr'),
            fitness_value=8,
        ),
    ),
)
def fitness_function(request):
    return request.param
