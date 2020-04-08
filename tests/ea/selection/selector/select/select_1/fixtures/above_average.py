import pytest
import stk

from .utilities import get_rank_fitness
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
population2 = (
    stk.MoleculeRecord(
        topology_graph=topology_graph,
    ).with_fitness_value(100),
    stk.MoleculeRecord(
        topology_graph=topology_graph,
    ).with_fitness_value(1),
)


@pytest.fixture(
    params=(
        CaseData(
            selector=stk.AboveAverage(),
            population=population1,
            selected=(
                stk.Batch(
                    records=(population1[0], ),
                    fitness_values={population1[0]: 10},
                    key_maker=stk.Inchi(),
                ),
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
            ),
        ),
        CaseData(
            selector=stk.AboveAverage(num_batches=2),
            population=population1,
            selected=(
                stk.Batch(
                    records=(population1[0], ),
                    fitness_values={population1[0]: 10},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[0], ),
                    fitness_values={population1[0]: 10},
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        CaseData(
            selector=stk.AboveAverage(
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
            ),
        ),
        CaseData(
            selector=stk.AboveAverage(
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
            ),
        ),
        CaseData(
            selector=stk.AboveAverage(
                batch_size=2,
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
                    records=(population1[0], population1[3], ),
                    fitness_values={
                        population1[0]: 10,
                        population1[3]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[0], population1[4], ),
                    fitness_values={
                        population1[0]: 10,
                        population1[4]: 1,
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
                    records=(population1[1], population1[3], ),
                    fitness_values={
                        population1[1]: 9,
                        population1[3]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[1], population1[4], ),
                    fitness_values={
                        population1[1]: 9,
                        population1[4]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        CaseData(
            selector=stk.AboveAverage(
                num_batches=3,
                batch_size=2,
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
            ),
        ),
        CaseData(
            selector=stk.AboveAverage(
                batch_size=2,
                duplicate_molecules=False,
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
            ),
        ),
        CaseData(
            selector=stk.AboveAverage(
                batch_size=2,
                duplicate_batches=False,
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
                    records=(population1[0], population1[3], ),
                    fitness_values={
                        population1[0]: 10,
                        population1[3]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[0], population1[4], ),
                    fitness_values={
                        population1[0]: 10,
                        population1[4]: 1,
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
                    records=(population1[1], population1[3], ),
                    fitness_values={
                        population1[1]: 9,
                        population1[3]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[1], population1[4], ),
                    fitness_values={
                        population1[1]: 9,
                        population1[4]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        CaseData(
            selector=stk.AboveAverage(
                fitness_modifier=get_rank_fitness,
            ),
            population=population2,
            selected=(
                stk.Batch(
                    records=(population2[0], ),
                    fitness_values={population2[0]: 1},
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
    ),
)
def above_average(request):
    return request.param
