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
            fitness_normalizer=stk.Multiply(
                coefficient=(1, 2, 3),
                filter=lambda population, record:
                    record.get_fitness_value() is not None,
            ),
            population=(
                stk.MoleculeRecord(
                    topology_graph=topology_graph,
                ).with_fitness_value((1, 10, 100)),
                stk.MoleculeRecord(
                    topology_graph=topology_graph,
                ).with_fitness_value((2, 20, 200)),
                stk.MoleculeRecord(
                    topology_graph=topology_graph,
                ).with_fitness_value((3, 30, 300)),
                stk.MoleculeRecord(topology_graph),
            ),
            normalized=(
                stk.MoleculeRecord(
                    topology_graph=topology_graph,
                ).with_fitness_value((1, 20, 300)),
                stk.MoleculeRecord(
                    topology_graph=topology_graph,
                ).with_fitness_value((2, 40, 600)),
                stk.MoleculeRecord(
                    topology_graph=topology_graph,
                ).with_fitness_value((3, 60, 900)),
                stk.MoleculeRecord(topology_graph),
            ),
        ),
    ),
)
def multiply(request):
    return request.param
