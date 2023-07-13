from typing import Any

import pandas as pd
import pytest
import stk

from .case_data import CaseData


def _get_topology_graph() -> stk.polymer.Linear:
    return stk.polymer.Linear(
        building_blocks=[
            stk.BuildingBlock("BrCCBr", stk.BromoFactory()),
        ],
        repeating_unit="A",
        num_repeating_units=2,
    )


def get_generation(
    *fitness_values: float,
) -> dict[stk.MoleculeRecord[Any], float | None]:
    v1, v2, v3, *_ = fitness_values
    topology_graph = _get_topology_graph()
    return {
        stk.MoleculeRecord(
            topology_graph=topology_graph,
        ): v1,
        stk.MoleculeRecord(
            topology_graph=topology_graph,
        ): v2,
        stk.MoleculeRecord(
            topology_graph=topology_graph,
        ): v3,
        stk.MoleculeRecord(
            topology_graph=topology_graph,
        ): None,
    }


@pytest.fixture(
    scope="session",
    params=(
        lambda: CaseData(
            plotter=stk.ProgressPlotter(
                property=[
                    (0, 1, 2),
                    (10, 20, 30),
                    (40, 50, 60),
                    (40, 50, 60),
                    (70, 80, 90),
                ],
                y_label="Fitness Value",
            ),
            plot_data=pd.DataFrame(
                {
                    "Generation": [0] * 3
                    + [1] * 3
                    + [2] * 3
                    + [3] * 3
                    + [4] * 3,
                    "Fitness Value": [
                        2.0,
                        1.0,
                        0.0,
                        30.0,
                        20.0,
                        10.0,
                        60.0,
                        50.0,
                        40.0,
                        60.0,
                        50.0,
                        40.0,
                        90.0,
                        80.0,
                        70.0,
                    ],
                    "Type": ["Max", "Mean", "Min"] * 5,
                }
            ),
        ),
    ),
)
def case_data(request: pytest.FixtureRequest) -> CaseData:
    return request.param()
