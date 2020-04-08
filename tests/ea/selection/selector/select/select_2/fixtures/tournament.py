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
    ).with_fitness_value(0),
)


@pytest.fixture(
    params=(
        CaseData(
            selector=stk.Tournament(
                duplicate_molecules=False,
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
        CaseData(
            selector=stk.Tournament(
                duplicate_batches=False,
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
        CaseData(
            selector=stk.Tournament(
                batch_size=2,
                duplicate_batches=False
            ),
            population=population1,
            selected=(
                stk.Batch(
                    records=(population1[0], population1[1]),
                    fitness_values={
                        population1[0]: 10,
                        population1[1]: 9,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[0], population1[2], ),
                    fitness_values={
                        population1[0]: 10,
                        population1[2]: 2,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[0], population1[3]),
                    fitness_values={
                        population1[0]: 10,
                        population1[3]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[0], population1[4]),
                    fitness_values={
                        population1[0]: 10,
                        population1[4]: 0,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[1], population1[2], ),
                    fitness_values={
                        population1[1]: 9,
                        population1[2]: 2,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[1], population1[3]),
                    fitness_values={
                        population1[1]: 9,
                        population1[3]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[1], population1[4]),
                    fitness_values={
                        population1[1]: 9,
                        population1[4]: 0,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[2], population1[3], ),
                    fitness_values={
                        population1[2]: 2,
                        population1[3]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[2], population1[4], ),
                    fitness_values={
                        population1[2]: 2,
                        population1[4]: 0,
                    },
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
    ),
)
def tournament(request):
    return request.param
