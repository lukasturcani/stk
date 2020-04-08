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
    ).with_fitness_value(11),
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
)
population2 = (
    stk.MoleculeRecord(
        topology_graph=get_topology_graph(6),
    ).with_fitness_value(100),
    stk.MoleculeRecord(
        topology_graph=get_topology_graph(7),
    ).with_fitness_value(100),
    stk.MoleculeRecord(
        topology_graph=get_topology_graph(8),
    ).with_fitness_value(1),
)


@pytest.fixture(
    params=(
        CaseData(
            selector=stk.Worst(),
            population=population1,
            selected=(
                stk.Batch(
                    records=(population1[4], ),
                    fitness_values={population1[4]: 1},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[3], ),
                    fitness_values={population1[3]: 2},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[2], ),
                    fitness_values={population1[2]: 9},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[1], ),
                    fitness_values={population1[1]: 10},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[0], ),
                    fitness_values={population1[0]: 11},
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        CaseData(
            selector=stk.Worst(num_batches=2),
            population=population1,
            selected=(
                stk.Batch(
                    records=(population1[4], ),
                    fitness_values={population1[4]: 1},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[3], ),
                    fitness_values={population1[3]: 2},
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        CaseData(
            selector=stk.Worst(
                duplicate_molecules=False,
            ),
            population=population2,
            selected=(
                stk.Batch(
                    records=(population2[2], ),
                    fitness_values={population2[2]: 1},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population2[0], ),
                    fitness_values={population2[0]: 100},
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        CaseData(
            selector=stk.Worst(
                duplicate_batches=False,
            ),
            population=population2,
            selected=(
                stk.Batch(
                    records=(population2[2], ),
                    fitness_values={population2[2]: 1},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population2[0], ),
                    fitness_values={population2[0]: 100},
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        CaseData(
            selector=stk.Worst(
                batch_size=2,
            ),
            population=population1,
            selected=(
                stk.Batch(
                    records=(population1[3], population1[4], ),
                    fitness_values={
                        population1[3]: 2,
                        population1[4]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[2], population1[4], ),
                    fitness_values={
                        population1[2]: 9,
                        population1[4]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[1], population1[4], ),
                    fitness_values={
                        population1[1]: 10,
                        population1[4]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[2], population1[3], ),
                    fitness_values={
                        population1[2]: 9,
                        population1[3]: 2,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[0], population1[4], ),
                    fitness_values={
                        population1[0]: 11,
                        population1[4]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[1], population1[3], ),
                    fitness_values={
                        population1[1]: 10,
                        population1[3]: 2,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[0], population1[3], ),
                    fitness_values={
                        population1[0]: 11,
                        population1[3]: 2,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[1], population1[2], ),
                    fitness_values={
                        population1[1]: 10,
                        population1[2]: 9,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[0], population1[2], ),
                    fitness_values={
                        population1[0]: 11,
                        population1[2]: 9,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[0], population1[1]),
                    fitness_values={
                        population1[0]: 11,
                        population1[1]: 10,
                    },
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        CaseData(
            selector=stk.Worst(
                num_batches=3,
                batch_size=2,
            ),
            population=population1,
            selected=(
                stk.Batch(
                    records=(population1[3], population1[4], ),
                    fitness_values={
                        population1[3]: 2,
                        population1[4]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[2], population1[4], ),
                    fitness_values={
                        population1[2]: 9,
                        population1[4]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[1], population1[4], ),
                    fitness_values={
                        population1[1]: 10,
                        population1[4]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        CaseData(
            selector=stk.Worst(
                batch_size=2,
                duplicate_molecules=False,
            ),
            population=population1,
            selected=(
                stk.Batch(
                    records=(population1[3], population1[4], ),
                    fitness_values={
                        population1[3]: 2,
                        population1[4]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[1], population1[2], ),
                    fitness_values={
                        population1[1]: 10,
                        population1[2]: 9,
                    },
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        CaseData(
            selector=stk.Worst(
                batch_size=2,
                duplicate_batches=False,
            ),
            population=population2,
            selected=(
                stk.Batch(
                    records=(population2[0], population2[2], ),
                    fitness_values={
                        population2[0]: 100,
                        population2[2]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population2[0], population2[1]),
                    fitness_values={
                        population2[0]: 100,
                        population2[1]: 100,
                    },
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
    ),
)
def worst(request):
    return request.param
