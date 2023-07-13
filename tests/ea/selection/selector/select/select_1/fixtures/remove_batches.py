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
def remove_batches_population_1() -> dict[stk.MoleculeRecord[Any], float]:
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
        ): 1,
    }


@pytest.fixture(
    scope="session",
    params=(
        lambda population: CaseData.new(
            selector=stk.RemoveBatches(
                remover=stk.Worst(4),
                selector=stk.Best(),
            ),
            population=population,
            selected=[[0]],
        ),
    ),
)
def remove_batches(
    request: pytest.FixtureRequest,
    remove_batches_population_1: dict[stk.MoleculeRecord[Any], float],
) -> CaseData:
    return request.param(remove_batches_population_1)
