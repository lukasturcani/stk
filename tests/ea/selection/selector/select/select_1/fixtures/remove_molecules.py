from __future__ import annotations

import pytest
import stk

from ..case_data import CaseData


def get_topology_graph(num_repeating_units):
    return stk.polymer.Linear(
        building_blocks=(stk.BuildingBlock("BrCCBr", [stk.BromoFactory()]),),
        repeating_unit="A",
        num_repeating_units=num_repeating_units,
    )


@pytest.fixture(scope="session")
def remove_molecules_population_1() -> tuple[stk.MoleculeRecord, ...]:
    return (
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
    scope="session",
    params=(
        lambda population: CaseData(
            selector=stk.RemoveMolecules(
                remover=stk.Best(2),
                selector=stk.Best(),
            ),
            population=population,
            selected=(
                stk.Batch(
                    records=(population[2],),
                    fitness_values={population[2]: 2},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population[3],),
                    fitness_values={population[3]: 1},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population[4],),
                    fitness_values={population[4]: 1},
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
    ),
)
def remove_molecules(
    request,
    remove_molecules_population_1: tuple[stk.MoleculeRecord, ...],
) -> CaseData:
    return request.param(remove_molecules_population_1)
