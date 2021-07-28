import pytest
import stk

from ..case_data import CaseData


def _get_case_data_1() -> CaseData:
    topology_graph = stk.polymer.Linear(
        building_blocks=(
            stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
        ),
        repeating_unit='A',
        num_repeating_units=2,
    )
    return CaseData(
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
                topology_graph=topology_graph,
            ).with_fitness_value(1),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ).with_fitness_value(2),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ).with_fitness_value(3),
            stk.MoleculeRecord(topology_graph),
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
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ).with_fitness_value(0.5),
        ),
    )


@pytest.fixture(
    scope='session',
    params=(
        _get_case_data_1,
    ),
)
def replace_fitness(request) -> CaseData:
    return request.param()
