import pytest
import stk

from ..case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            fitness_normalizer=stk.DivideByMean(
                filter=lambda population, record:
                    record.get_fitness_value() is not None,
            ),
            population=(
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCBr'),
                ).with_fitness_value((1, 10, 100)),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCCBr'),
                ).with_fitness_value((2, 20, 200)),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCCCBr'),
                ).with_fitness_value((3, 30, 300)),
                stk.MoleculeRecord(stk.BuildingBlock('BrCCCCBr')),
            ),
            normalized=(
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCBr'),
                ).with_fitness_value((0.5, 0.5, 0.5)),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCCBr'),
                ).with_fitness_value((1, 1, 1)),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCCCBr'),
                ).with_fitness_value((1.5, 1.5, 1.5)),
                stk.MoleculeRecord(stk.BuildingBlock('BrCCCCBr')),
            ),
        ),
    ),
)
def multiply(request):
    return request.param
