import pytest
import stk

from ..case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            fitness_normalizer=stk.Sum(
                filter=lambda population, record:
                    record.get_fitness_value() is not None,
            ),
            population=(
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCBr'),
                ).with_fitness_value((1, -5, 5)),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCCBr'),
                ).with_fitness_value((3, -10, 2)),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCCCBr'),
                ).with_fitness_value((2, 20, 1)),
                stk.MoleculeRecord(stk.BuildingBlock('BrCCCCBr')),
            ),
            normalized=(
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCBr'),
                ).with_fitness_value(1),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCCBr'),
                ).with_fitness_value(-5),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCCCBr'),
                ).with_fitness_value(23),
                stk.MoleculeRecord(stk.BuildingBlock('BrCCCCBr')),
            ),
        ),
    ),
)
def shift_up(request):
    return request.param
