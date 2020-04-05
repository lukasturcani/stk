import pytest
import stk

from ..case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            selector=stk.AboveAverage(),
            population=(
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCBr'),
                ).with_fitness_value(5),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCBr'),
                ).with_fitness_value(1),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCCBr'),
                ).with_fitness_value(1),
            ),
            selected=(
                stk.Batch(
                    records=(
                    ),
                    fitness_values={},
                    key_maker=stk.Inchi(),
                )
            ),
        ),
    ),
)
def above_average(request):
    return request.param
