import pytest
import stk

from ..case_data import CaseData


def _get_case_data_1() -> CaseData:
    topology_graph = stk.polymer.Linear(
        building_blocks=(
            stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
        ),
        repeating_unit='A',
        num_repeating_units=2,
    )
    return CaseData(
        fitness_normalizer=stk.ShiftUp(
            filter=lambda population, record:
                record.get_fitness_value() is not None,
        ),
        population=(
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ).with_fitness_value((1, -5, 5)),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ).with_fitness_value((3, -10, 2)),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ).with_fitness_value((2, 20, 1)),
            stk.MoleculeRecord(topology_graph),
        ),
        normalized=(
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ).with_fitness_value((1, 6, 5)),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ).with_fitness_value((3, 1, 2)),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ).with_fitness_value((2, 31, 1)),
            stk.MoleculeRecord(topology_graph),
        ),
    )


@pytest.fixture(
    scope='session',
    params=(
        _get_case_data_1,
    ),
)
def shift_up(request):
    return request.param()
