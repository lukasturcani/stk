import pytest
import stk

from ..case_data import CaseData


@pytest.fixture(
    scope="session",
    params=(
        lambda: CaseData(
            fitness_calculator=stk.FitnessFunction(
                fitness_function=stk.Molecule.get_num_atoms,
            ),
            molecule=stk.BuildingBlock("BrCCBr"),
            fitness_value=8,
        ),
    ),
)
def fitness_function(request) -> CaseData:
    return request.param()
