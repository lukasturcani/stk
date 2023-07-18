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
        fitness_normalizer=stk.ReplaceFitness(
            get_replacement=lambda fitness_values: min(
                fitness_value
                for fitness_value in fitness_values.values()
                if fitness_value is not None
            )
            / 2,
            filter=lambda fitness_values, record: fitness_values[record]
            is None,
        ),
        fitness_values={
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): (1, 1),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): (2, 2),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): (3, 3),
            stk.MoleculeRecord(topology_graph): (None, 0.5),
        },
    )


@pytest.fixture(
    scope="session",
    params=(_get_case_data_1,),
)
def replace_fitness(request: pytest.FixtureRequest) -> CaseData:
    return request.param()
