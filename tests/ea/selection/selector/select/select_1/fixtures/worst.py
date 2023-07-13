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
def worst_population_1() -> dict[stk.MoleculeRecord[Any], float]:
    return {
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(2),
        ): 11,
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(3),
        ): 10,
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(4),
        ): 9,
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(5),
        ): 2,
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(6),
        ): 1,
    }


@pytest.fixture(scope="session")
def worst_population_2() -> dict[stk.MoleculeRecord[Any], float]:
    return {
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(6),
        ): 100,
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(6),
        ): 100,
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(8),
        ): 1,
    }


@pytest.fixture(
    scope="session",
    params=(
        lambda population: CaseData.new(
            selector=stk.Worst(),
            population=population,
            selected=([4], [3], [2], [1], [0]),
        ),
        lambda population: CaseData.new(
            selector=stk.Worst(num_batches=2),
            population=population,
            selected=([4], [3]),
        ),
        lambda population: CaseData.new(
            selector=stk.Worst(
                batch_size=2,
            ),
            population=population,
            selected=(
                [3, 4],
                [2, 4],
                [1, 4],
                [2, 3],
                [0, 4],
                [1, 3],
                [0, 3],
                [1, 2],
                [0, 2],
                [0, 1],
            ),
        ),
        lambda population: CaseData.new(
            selector=stk.Worst(
                num_batches=3,
                batch_size=2,
            ),
            population=population,
            selected=([3, 4], [2, 4], [1, 4]),
        ),
        lambda population: CaseData.new(
            selector=stk.Worst(
                batch_size=2,
                duplicate_molecules=False,
            ),
            population=population,
            selected=([3, 4], [1, 2]),
        ),
    ),
)
def worst_population_1_case_data(
    request: pytest.FixtureRequest,
    worst_population_1: dict[stk.MoleculeRecord[Any], float],
) -> CaseData:
    return request.param(worst_population_1)


@pytest.fixture(
    scope="session",
    params=(
        lambda population: CaseData.new(
            selector=stk.Worst(
                duplicate_molecules=False,
            ),
            population=population,
            selected=([2], [0]),
        ),
        lambda population: CaseData.new(
            selector=stk.Worst(
                duplicate_batches=False,
            ),
            population=population,
            selected=([2], [0]),
        ),
        lambda population: CaseData.new(
            selector=stk.Worst(
                batch_size=2,
                duplicate_batches=False,
            ),
            population=population,
            selected=([0, 2], [0, 1]),
        ),
    ),
)
def worst_population_2_case_data(
    request: pytest.FixtureRequest,
    worst_population_2: dict[stk.MoleculeRecord[Any], float],
) -> CaseData:
    return request.param(worst_population_2)


@pytest.fixture(
    scope="session",
    params=(
        lazy_fixture("worst_population_1_case_data"),
        lazy_fixture("worst_population_2_case_data"),
    ),
)
def worst(request: pytest.FixtureRequest) -> CaseData:
    return request.param
