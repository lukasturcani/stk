import pytest
import stk

from ..case_data import CaseData

topology_graph = stk.polymer.Linear(
    building_blocks=(
        stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
    ),
    repeating_unit='A',
    num_repeating_units=2,
)


@pytest.fixture(
    params=(
        CaseData(
            fitness_normalizer=stk.NullFitnessNormalizer(),
            population=(
                stk.MoleculeRecord(
                    topology_graph=topology_graph,
                ).with_fitness_value(1),
                stk.MoleculeRecord(
                    topology_graph=topology_graph,
                ).with_fitness_value(2),
                stk.MoleculeRecord(
                    topology_graph=topology_graph,
                ).with_fitness_value(3),
            ),
            normalized=(
                stk.MoleculeRecord(
                    topology_graph=topology_graph,
                ).with_fitness_value(1),
                stk.MoleculeRecord(
                    topology_graph=topology_graph,
                ).with_fitness_value(2),
                stk.MoleculeRecord(
                    topology_graph=topology_graph,
                ).with_fitness_value(3),
            ),
        ),
    ),
)
def null(request):
    return request.param
