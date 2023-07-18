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
        fitness_normalizer=stk.NormalizerSequence(
            fitness_normalizers=(
                stk.Multiply(
                    coefficient=(1, 2, 4),
                    filter=lambda fitness_values, record: fitness_values[
                        record
                    ]
                    is not None,
                ),
                stk.Sum(
                    filter=lambda fitness_values, record: fitness_values[
                        record
                    ]
                    is not None,
                ),
            ),
        ),
        fitness_values={
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): ((4, 10, 1), 28),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): ((9, 20, 2), 57),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): ((16, 30, 4), 92),
            stk.MoleculeRecord(topology_graph): (None, None),
        },
    )


@pytest.fixture(
    scope="session",
    params=(_get_case_data_1,),
)
def sequence(request: pytest.FixtureRequest) -> CaseData:
    return request.param()
