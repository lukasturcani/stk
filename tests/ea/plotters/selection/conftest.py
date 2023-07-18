from typing import Any

import pytest
import stk

from .case_data import CaseData


def get_topology_graph(num_repeating_units: int) -> stk.TopologyGraph:
    return stk.polymer.Linear(
        building_blocks=[
            stk.BuildingBlock("BrCCBr", stk.BromoFactory()),
        ],
        repeating_unit="A",
        num_repeating_units=num_repeating_units,
    )


@pytest.fixture(
    scope="session",
    params=(
        lambda: CaseData(  # type: ignore
            selector=stk.Roulette(num_batches=50),
            population={
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(2),
                ): 1,
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(3),
                ): 2,
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(4),
                ): 3,
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(5),
                ): 4,
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(6),
                ): 5,
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(7),
                ): 6,
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(8),
                ): 7,
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(9),
                ): 8,
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(10),
                ): 9,
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(11),
                ): 10,
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(12),
                ): 11,
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(13),
                ): 12,
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(14),
                ): 13,
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(15),
                ): 14,
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(16),
                ): 15,
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(17),
                ): 16,
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(18),
                ): 17,
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(19),
                ): 18,
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(20),
                ): 19,
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(21),
                ): 20,
            },
        ),
    ),
)
def case_data(
    request: pytest.FixtureRequest,
) -> CaseData[stk.MoleculeRecord[Any]]:
    return request.param()
