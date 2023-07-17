from typing import Any

import pytest
import stk
from pytest_lazyfixture import lazy_fixture

from ..case_data import CaseData


def get_topology_graph(num_repeating_units: int) -> stk.TopologyGraph:
    return stk.polymer.Linear(
        building_blocks=[
            stk.BuildingBlock("BrCCBr", stk.BromoFactory()),
        ],
        repeating_unit="A",
        num_repeating_units=num_repeating_units,
    )


@pytest.fixture(scope="session")
def best_population_1() -> dict[stk.MoleculeRecord[Any], float]:
    return {
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(1),
        ): 1,
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(2),
        ): 10,
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(3),
        ): 9,
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(4),
        ): 4,
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(5),
        ): 3,
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(6),
        ): 2,
    }


@pytest.fixture(scope="session")
def best_population_2() -> dict[stk.MoleculeRecord[Any], float]:
    return {
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(7),
        ): 100,
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(7),
        ): 100,
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(9),
        ): 1,
    }


@pytest.fixture(
    scope="session",
    params=(
        lambda population: CaseData.new(
            selector=stk.Best(),
            population=population,
            selected=([1], [2], [3], [4], [5], [0]),
        ),
        lambda population: CaseData.new(
            selector=stk.Best(num_batches=2),
            population=population,
            selected=([1], [2]),
        ),
        lambda population: CaseData.new(
            selector=stk.Best(
                batch_size=2,
            ),
            population=population,
            selected=(
                [1, 2],
                [1, 3],
                [1, 4],
                [2, 3],
                [1, 5],
                [2, 4],
                [0, 1],
                [2, 5],
                [0, 2],
                [3, 4],
                [3, 5],
                [0, 3],
                [4, 5],
                [0, 4],
                [0, 5],
            ),
        ),
        lambda population: CaseData.new(
            selector=stk.Best(
                num_batches=3,
                batch_size=2,
            ),
            population=population,
            selected=([1, 2], [1, 3], [1, 4]),
        ),
        lambda population: CaseData.new(
            selector=stk.Best(
                batch_size=2,
                duplicate_molecules=False,
            ),
            population=population,
            selected=([1, 2], [3, 4], [0, 5]),
        ),
    ),
)
def best_population_1_case_data(
    request: pytest.FixtureRequest,
    best_population_1: dict[stk.MoleculeRecord[Any], float],
) -> CaseData:
    return request.param(best_population_1)


@pytest.fixture(
    scope="session",
    params=(
        lambda population: CaseData.new(
            selector=stk.Best(
                duplicate_molecules=False,
            ),
            population=population,
            selected=([0], [2]),
        ),
        lambda population: CaseData.new(
            selector=stk.Best(
                duplicate_batches=False,
            ),
            population=population,
            selected=([0], [2]),
        ),
        lambda population: CaseData.new(
            selector=stk.Best(
                batch_size=2,
                duplicate_batches=False,
            ),
            population=population,
            selected=([0, 1], [0, 2]),
        ),
    ),
)
def best_population_2_case_data(
    request: pytest.FixtureRequest,
    best_population_2: dict[stk.MoleculeRecord[Any], float],
) -> CaseData:
    return request.param(best_population_2)


@pytest.fixture(
    scope="session",
    params=(
        lazy_fixture("best_population_1_case_data"),
        lazy_fixture("best_population_2_case_data"),
    ),
)
def best(request: pytest.FixtureRequest) -> CaseData:
    return request.param
