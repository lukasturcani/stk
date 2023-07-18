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
        fitness_normalizer=stk.Power(
            power=(0.5, 0, -1),
            filter=lambda fitness_values, record: fitness_values[record]
            is not None,
        ),
        fitness_values={
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): ((4, 10, 1), (2, 1, 1)),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): ((9, 20, 2), (3, 1, 0.5)),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): ((16, 30, 4), (4, 1, 0.25)),
            stk.MoleculeRecord(topology_graph): (None, None),
        },
    )


@pytest.fixture(
    scope="session",
    params=(_get_case_data_1,),
)
def power(request: pytest.FixtureRequest) -> CaseData:
    return request.param()
