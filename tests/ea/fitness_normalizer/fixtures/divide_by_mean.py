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
            fitness_normalizer=stk.DivideByMean(),
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
                ).with_fitness_value(0.5),
                stk.MoleculeRecord(
                    topology_graph=topology_graph,
                ).with_fitness_value(1),
                stk.MoleculeRecord(
                    topology_graph=topology_graph,
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
                ).with_fitness_value((0.5, 0.5, 0.5)),
                stk.MoleculeRecord(
                    topology_graph=topology_graph,
                ).with_fitness_value((1, 1, 1)),
                stk.MoleculeRecord(
                    topology_graph=topology_graph,
                ).with_fitness_value((1.5, 1.5, 1.5)),
                stk.MoleculeRecord(topology_graph),
            ),
        ),
    ),
)
def divide_by_mean(request):
    return request.param
