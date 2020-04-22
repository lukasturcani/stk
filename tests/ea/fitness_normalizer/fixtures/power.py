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
            fitness_normalizer=stk.Power(
                power=(0.5, 0, -1),
                filter=lambda population, record:
                    record.get_fitness_value() is not None,
            ),
            population=(
                stk.MoleculeRecord(
                    topology_graph=topology_graph,
                ).with_fitness_value((4, 10, 1)),
                stk.MoleculeRecord(
                    topology_graph=topology_graph,
                ).with_fitness_value((9, 20, 2)),
                stk.MoleculeRecord(
                    topology_graph=topology_graph,
                ).with_fitness_value((16, 30, 4)),
                stk.MoleculeRecord(topology_graph),
            ),
            normalized=(
                stk.MoleculeRecord(
                    topology_graph=topology_graph,
                ).with_fitness_value((2, 1, 1)),
                stk.MoleculeRecord(
                    topology_graph=topology_graph,
                ).with_fitness_value((3, 1, 0.5)),
                stk.MoleculeRecord(
                    topology_graph=topology_graph,
                ).with_fitness_value((4, 1, 0.25)),
                stk.MoleculeRecord(topology_graph),
            ),
        ),
    ),
)
def power(request):
    return request.param
