import pytest
import stk

from ..case_data import CaseData


population1 = (
    stk.MoleculeRecord(
        molecule=stk.BuildingBlock('BrCBr'),
    ).with_fitness_value(5),
    stk.MoleculeRecord(
        molecule=stk.BuildingBlock('BrCCBr'),
    ).with_fitness_value(1),
    stk.MoleculeRecord(
        molecule=stk.BuildingBlock('BrCCCBr'),
    ).with_fitness_value(1),
)


@pytest.fixture(
    params=(
        CaseData(
            selector=stk.AboveAverage(),
            population=population1,
            selected=(
                stk.Batch(
                    records=(
                        population1[0],
                    ),
                    fitness_values={
                        population1[0]: 5,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(
                        population1[0]
                    ),
                    fitness_values={
                        population1[0]: 5,
                    },
                    key_maker=stk.Inchi(),
                )
            ),
        ),
    ),
)
def above_average(request):
    return request.param
