from __future__ import annotations

import pytest
from pytest_lazyfixture import lazy_fixture
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


@pytest.fixture(scope='session')
def worst_population_1() -> tuple[stk.MoleculeRecord, ...]:
    return (
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(2),
        ).with_fitness_value(11),
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(3),
        ).with_fitness_value(10),
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(4),
        ).with_fitness_value(9),
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(5),
        ).with_fitness_value(2),
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(6),
        ).with_fitness_value(1),
    )


@pytest.fixture(scope='session')
def worst_population_2() -> tuple[stk.MoleculeRecord, ...]:
    return (
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(6),
        ).with_fitness_value(100),
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(6),
        ).with_fitness_value(100),
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(8),
        ).with_fitness_value(1),
    )


@pytest.fixture(
    scope='session',
    params=(
        lambda population: CaseData(
            selector=stk.Worst(),
            population=population,
            selected=(
                stk.Batch(
                    records=(population[4], ),
                    fitness_values={population[4]: 1},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population[3], ),
                    fitness_values={population[3]: 2},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population[2], ),
                    fitness_values={population[2]: 9},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population[1], ),
                    fitness_values={population[1]: 10},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population[0], ),
                    fitness_values={population[0]: 11},
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        lambda population: CaseData(
            selector=stk.Worst(num_batches=2),
            population=population,
            selected=(
                stk.Batch(
                    records=(population[4], ),
                    fitness_values={population[4]: 1},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population[3], ),
                    fitness_values={population[3]: 2},
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        lambda population: CaseData(
            selector=stk.Worst(
                batch_size=2,
            ),
            population=population,
            selected=(
                stk.Batch(
                    records=(population[3], population[4], ),
                    fitness_values={
                        population[3]: 2,
                        population[4]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population[2], population[4], ),
                    fitness_values={
                        population[2]: 9,
                        population[4]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population[1], population[4], ),
                    fitness_values={
                        population[1]: 10,
                        population[4]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population[2], population[3], ),
                    fitness_values={
                        population[2]: 9,
                        population[3]: 2,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population[0], population[4], ),
                    fitness_values={
                        population[0]: 11,
                        population[4]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population[1], population[3], ),
                    fitness_values={
                        population[1]: 10,
                        population[3]: 2,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population[0], population[3], ),
                    fitness_values={
                        population[0]: 11,
                        population[3]: 2,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population[1], population[2], ),
                    fitness_values={
                        population[1]: 10,
                        population[2]: 9,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population[0], population[2], ),
                    fitness_values={
                        population[0]: 11,
                        population[2]: 9,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population[0], population[1]),
                    fitness_values={
                        population[0]: 11,
                        population[1]: 10,
                    },
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        lambda population: CaseData(
            selector=stk.Worst(
                num_batches=3,
                batch_size=2,
            ),
            population=population,
            selected=(
                stk.Batch(
                    records=(population[3], population[4], ),
                    fitness_values={
                        population[3]: 2,
                        population[4]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population[2], population[4], ),
                    fitness_values={
                        population[2]: 9,
                        population[4]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population[1], population[4], ),
                    fitness_values={
                        population[1]: 10,
                        population[4]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        lambda population: CaseData(
            selector=stk.Worst(
                batch_size=2,
                duplicate_molecules=False,
            ),
            population=population,
            selected=(
                stk.Batch(
                    records=(population[3], population[4], ),
                    fitness_values={
                        population[3]: 2,
                        population[4]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population[1], population[2], ),
                    fitness_values={
                        population[1]: 10,
                        population[2]: 9,
                    },
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
    ),
)
def worst_population_1_case_data(
    request,
    worst_population_1: tuple[stk.MoleculeRecord, ...],
) -> CaseData:
    return request.param(worst_population_1)


@pytest.fixture(
    scope='session',
    params=(
        lambda population: CaseData(
            selector=stk.Worst(
                duplicate_molecules=False,
            ),
            population=population,
            selected=(
                stk.Batch(
                    records=(population[2], ),
                    fitness_values={population[2]: 1},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population[0], ),
                    fitness_values={population[0]: 100},
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        lambda population: CaseData(
            selector=stk.Worst(
                duplicate_batches=False,
            ),
            population=population,
            selected=(
                stk.Batch(
                    records=(population[2], ),
                    fitness_values={population[2]: 1},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population[0], ),
                    fitness_values={population[0]: 100},
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        lambda population: CaseData(
            selector=stk.Worst(
                batch_size=2,
                duplicate_batches=False,
            ),
            population=population,
            selected=(
                stk.Batch(
                    records=(population[0], population[2], ),
                    fitness_values={
                        population[0]: 100,
                        population[2]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population[0], population[1]),
                    fitness_values={
                        population[0]: 100,
                        population[1]: 100,
                    },
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
    ),
)
def worst_population_2_case_data(
    request,
    worst_population_2: tuple[stk.MoleculeRecord, ...],
) -> CaseData:
    return request.param(worst_population_2)


@pytest.fixture(
    scope='session',
    params=(
        lazy_fixture('worst_population_1_case_data'),
        lazy_fixture('worst_population_2_case_data'),
    ),
)
def worst(request):
    return request.param
