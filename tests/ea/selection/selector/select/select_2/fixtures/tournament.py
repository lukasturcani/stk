from typing import Any

import pytest
import stk

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
def tournament_population_1() -> dict[stk.MoleculeRecord[Any], float]:
    return {
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(2),
        ): 10,
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(3),
        ): 9,
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(4),
        ): 2,
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(5),
        ): 1,
        stk.MoleculeRecord(
            topology_graph=get_topology_graph(6),
        ): 0,
    }


@pytest.fixture(
    scope="session",
    params=(
        lambda population: CaseData.new(
            selector=stk.Tournament(
                duplicate_molecules=False,
            ),
            population=population,
            selected=([0], [1], [2], [3]),
        ),
        lambda population: CaseData.new(
            selector=stk.Tournament(
                duplicate_batches=False,
            ),
            population=population,
            selected=([0], [1], [2], [3]),
        ),
        lambda population: CaseData.new(
            selector=stk.Tournament(batch_size=2, duplicate_batches=False),
            population=population,
            selected=(
                [0, 1],
                [0, 2],
                [0, 3],
                [0, 4],
                [1, 2],
                [1, 3],
                [1, 4],
                [2, 3],
                [2, 4],
            ),
        ),
    ),
)
def tournament(
    request: pytest.FixtureRequest,
    tournament_population_1: dict[stk.MoleculeRecord[Any], float],
) -> CaseData:
    return request.param(tournament_population_1)
