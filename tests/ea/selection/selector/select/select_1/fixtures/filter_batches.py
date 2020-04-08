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


population1 = (
    stk.MoleculeRecord(
        topology_graph=topology_graph,
    ).with_fitness_value(10),
    stk.MoleculeRecord(
        topology_graph=topology_graph,
    ).with_fitness_value(9),
    stk.MoleculeRecord(
        topology_graph=topology_graph,
    ).with_fitness_value(2),
    stk.MoleculeRecord(
        topology_graph=topology_graph,
    ).with_fitness_value(1),
    stk.MoleculeRecord(
        topology_graph=topology_graph,
    ).with_fitness_value(1),
)


@pytest.fixture(
    params=(
        CaseData(
            selector=stk.FilterBatches(
                filter=stk.Best(4),
                selector=stk.Best(),
            ),
            population=population1,
            selected=(
                stk.Batch(
                    records=(population1[0], ),
                    fitness_values={population1[0]: 10},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[1], ),
                    fitness_values={population1[1]: 9},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[2], ),
                    fitness_values={population1[2]: 2},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[3], ),
                    fitness_values={population1[3]: 1},
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
    ),
)
def filter_batches(request):
    return request.param
