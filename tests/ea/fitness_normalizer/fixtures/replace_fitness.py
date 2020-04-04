import pytest
import stk

from ..case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            fitness_normalizer=stk.ReplaceFitness(
                get_replacement=lambda population:
                    min(
                        record.get_fitness_value()
                        for record in population
                        if record.get_fitness_value() is not None
                    ),
                filter=lambda population, record:
                    record.get_fitness_value() is None,
            ),
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
                stk.MoleculeRecord(stk.BuildingBlock('BrCCCCBr')),
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
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCCCBr'),
                ).with_fitness_value(0.5),
            ),
        ),
    ),
)
def replace_fitness(request):
    return request.param
