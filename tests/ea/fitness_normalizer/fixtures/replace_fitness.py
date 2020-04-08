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
                    )/2,
                filter=lambda population, record:
                    record.get_fitness_value() is None,
            ),
            population=(
                stk.MoleculeRecord(
                    graph=None,
                ).with_fitness_value(1),
                stk.MoleculeRecord(
                    graph=None,
                ).with_fitness_value(2),
                stk.MoleculeRecord(
                    graph=None,
                ).with_fitness_value(3),
                stk.MoleculeRecord(None),
            ),
            normalized=(
                stk.MoleculeRecord(
                    graph=None,
                ).with_fitness_value(1),
                stk.MoleculeRecord(
                    graph=None,
                ).with_fitness_value(2),
                stk.MoleculeRecord(
                    graph=None,
                ).with_fitness_value(3),
                stk.MoleculeRecord(
                    graph=None,
                ).with_fitness_value(0.5),
            ),
        ),
    ),
)
def replace_fitness(request):
    return request.param
