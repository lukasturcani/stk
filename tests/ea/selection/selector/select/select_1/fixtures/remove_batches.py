import pytest
import stk

from ..case_data import CaseData


def get_topology_graph(num_repeating_units):
    return stk.polymer.Linear(
        building_blocks=(
            stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
        ),
        repeating_unit='A',
        num_repeating_units=num_repeating_units,
    )


population1 = (
    stk.MoleculeRecord(
        topology_graph=get_topology_graph(2),
    ).with_fitness_value(10),
    stk.MoleculeRecord(
        topology_graph=get_topology_graph(3),
    ).with_fitness_value(9),
    stk.MoleculeRecord(
        topology_graph=get_topology_graph(4),
    ).with_fitness_value(2),
    stk.MoleculeRecord(
        topology_graph=get_topology_graph(5),
    ).with_fitness_value(1),
    stk.MoleculeRecord(
        topology_graph=get_topology_graph(6),
    ).with_fitness_value(1),
)


@pytest.fixture(
    params=(
        CaseData(
            selector=stk.RemoveBatches(
                remover=stk.Worst(4),
                selector=stk.Best(),
            ),
            population=population1,
            selected=(
                stk.Batch(
                    records=(population1[0], ),
                    fitness_values={population1[0]: 10},
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
    ),
)
def remove_batches(request):
    return request.param
