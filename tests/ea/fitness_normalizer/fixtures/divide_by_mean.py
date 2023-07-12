import pytest
import stk

from ..case_data import CaseData


def _get_case_data_1() -> CaseData:
    topology_graph = stk.polymer.Linear(
        building_blocks=[
            stk.BuildingBlock("BrCCBr", stk.BromoFactory()),
        ],
        repeating_unit="A",
        num_repeating_units=2,
    )
    return CaseData.new(
        fitness_normalizer=stk.DivideByMean(),
        fitness_values={
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): (1, 0.5),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): (2, 1),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): (3, 1.5),
        },
    )


def _get_case_data_2() -> CaseData:
    topology_graph = stk.polymer.Linear(
        building_blocks=[
            stk.BuildingBlock("BrCCBr", stk.BromoFactory()),
        ],
        repeating_unit="A",
        num_repeating_units=2,
    )
    return CaseData.new(
        fitness_normalizer=stk.DivideByMean(
            filter=lambda fitness_values, record: fitness_values[record]
            is not None,
        ),
        fitness_values={
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): ((1, 10, 100), (0.5, 0.5, 0.5)),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): ((2, 20, 200), (1, 1, 1)),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): ((3, 30, 300), (1.5, 1.5, 1.5)),
            stk.MoleculeRecord(topology_graph): (None, None),
        },
    )


@pytest.fixture(
    scope="session",
    params=(
        _get_case_data_1,
        _get_case_data_2,
    ),
)
def divide_by_mean(request: pytest.FixtureRequest) -> CaseData:
    return request.param()
