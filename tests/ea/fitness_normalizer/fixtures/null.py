import pytest
import stk

from ..case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            fitness_normalizer=stk.Null(),
            population=(
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCBr'),
                ).with_fitness_value(1),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCCBr'),
                ).with_fitness_value(2),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCCCBr'),
                ).with_fitness_value(3),
            ),
            normalized=(
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCBr'),
                ).with_fitness_value(1),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCCBr'),
                ).with_fitness_value(2),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCCCBr'),
                ).with_fitness_value(3),
            ),
        ),
    ),
)
def null(request):
    return request.param
