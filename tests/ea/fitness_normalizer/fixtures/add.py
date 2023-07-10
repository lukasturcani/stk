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
        fitness_normalizer=stk.Add(
            number=(1, 2, 3),
            filter=lambda population, record: record.get_fitness_value()
            is not None,
        ),
        population={
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): (1, 10, 100),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): (2, 20, 200),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): (3, 30, 300),
            stk.MoleculeRecord(topology_graph): None,
        },
        normalized={
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): (2, 12, 103),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): (3, 22, 203),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ): (4, 32, 303),
            stk.MoleculeRecord(topology_graph): None,
        },
    )


@pytest.fixture(
    scope="session",
    params=(_get_case_data_1,),
)
def add(request: pytest.FixtureRequest) -> CaseData:
    return request.param()
