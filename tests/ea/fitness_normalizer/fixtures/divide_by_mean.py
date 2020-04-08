import pytest
import stk

from ..case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            fitness_normalizer=stk.DivideByMean(),
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
            ),
            normalized=(
                stk.MoleculeRecord(
                    graph=None,
                ).with_fitness_value(0.5),
                stk.MoleculeRecord(
                    graph=None,
                ).with_fitness_value(1),
                stk.MoleculeRecord(
                    graph=None,
                ).with_fitness_value(1.5),
            ),
        ),
        CaseData(
            fitness_normalizer=stk.DivideByMean(
                filter=lambda population, record:
                    record.get_fitness_value() is not None,
            ),
            population=(
                stk.MoleculeRecord(
                    graph=None,
                ).with_fitness_value((1, 10, 100)),
                stk.MoleculeRecord(
                    graph=None,
                ).with_fitness_value((2, 20, 200)),
                stk.MoleculeRecord(
                    graph=None,
                ).with_fitness_value((3, 30, 300)),
                stk.MoleculeRecord(None),
            ),
            normalized=(
                stk.MoleculeRecord(
                    graph=None,
                ).with_fitness_value((0.5, 0.5, 0.5)),
                stk.MoleculeRecord(
                    graph=None,
                ).with_fitness_value((1, 1, 1)),
                stk.MoleculeRecord(
                    graph=None,
                ).with_fitness_value((1.5, 1.5, 1.5)),
                stk.MoleculeRecord(None),
            ),
        ),
    ),
)
def divide_by_mean(request):
    return request.param
