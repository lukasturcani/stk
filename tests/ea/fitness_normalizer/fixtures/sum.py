import pytest
import stk

from ..case_data import CaseData


def _get_case_data_1() -> CaseData:
    topology_graph = stk.polymer.Linear(
        building_blocks=(stk.BuildingBlock("BrCCBr", [stk.BromoFactory()]),),
        repeating_unit="A",
        num_repeating_units=2,
    )
    return CaseData(
        fitness_normalizer=stk.Sum(
            filter=lambda fitness_values, record: fitness_values[record]
            is not None,
        ),
        fitness_values={
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): (1, -5, 5),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): (3, -10, 2),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): (2, 20, 1),
            stk.MoleculeRecord(topology_graph): None,
        },
        normalized={
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): 1,
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): -5,
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): 23,
            stk.MoleculeRecord(topology_graph): None,
        },
    )


@pytest.fixture(
    scope="session",
    params=(_get_case_data_1,),
)
def sum(request: pytest.FixtureRequest) -> CaseData:
    return request.param()
