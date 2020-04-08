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
                    graph=None,
                ).with_fitness_value((1, -5, 5)),
                stk.MoleculeRecord(
                    graph=None,
                ).with_fitness_value((3, -10, 2)),
                stk.MoleculeRecord(
                    graph=None,
                ).with_fitness_value((2, 20, 1)),
                stk.MoleculeRecord(None),
            ),
            normalized=(
                stk.MoleculeRecord(
                    graph=None,
                ).with_fitness_value(1),
                stk.MoleculeRecord(
                    graph=None,
                ).with_fitness_value(-5),
                stk.MoleculeRecord(
                    graph=None,
                ).with_fitness_value(23),
                stk.MoleculeRecord(None),
            ),
        ),
    ),
)
def sum(request):
    return request.param
