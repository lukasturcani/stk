import pytest
import stk

from ..case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            fitness_normalizer=stk.Multiply(
                coefficient=(1, 2, 3),
                filter=lambda population, record:
                    record.get_fitness_value() is not None,
            ),
            population=(
                stk.MoleculeRecord(
                    topology_graph=None,
                ).with_fitness_value((1, 10, 100)),
                stk.MoleculeRecord(
                    topology_graph=None,
                ).with_fitness_value((2, 20, 200)),
                stk.MoleculeRecord(
                    topology_graph=None,
                ).with_fitness_value((3, 30, 300)),
                stk.MoleculeRecord(None),
            ),
            normalized=(
                stk.MoleculeRecord(
                    topology_graph=None,
                ).with_fitness_value((1, 20, 300)),
                stk.MoleculeRecord(
                    topology_graph=None,
                ).with_fitness_value((2, 40, 600)),
                stk.MoleculeRecord(
                    topology_graph=None,
                ).with_fitness_value((3, 60, 900)),
                stk.MoleculeRecord(None),
            ),
        ),
    ),
)
def multiply(request):
    return request.param
