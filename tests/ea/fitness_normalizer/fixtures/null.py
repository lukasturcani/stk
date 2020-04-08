import pytest
import stk

from ..case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            fitness_normalizer=stk.NullFitnessNormalizer(),
            population=(
                stk.MoleculeRecord(
                    topology_graph=None,
                ).with_fitness_value(1),
                stk.MoleculeRecord(
                    topology_graph=None,
                ).with_fitness_value(2),
                stk.MoleculeRecord(
                    topology_graph=None,
                ).with_fitness_value(3),
            ),
            normalized=(
                stk.MoleculeRecord(
                    topology_graph=None,
                ).with_fitness_value(1),
                stk.MoleculeRecord(
                    topology_graph=None,
                ).with_fitness_value(2),
                stk.MoleculeRecord(
                    topology_graph=None,
                ).with_fitness_value(3),
            ),
        ),
    ),
)
def null(request):
    return request.param
