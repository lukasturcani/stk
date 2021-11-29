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
        fitness_normalizer=stk.NormalizerSequence(
            fitness_normalizers=(
                stk.Multiply(
                    coefficient=(1, 2, 4),
                    filter=lambda population, record:
                        record.get_fitness_value() is not None,
                ),
                stk.Sum(
                    filter=lambda population, record:
                        record.get_fitness_value() is not None,
                ),
            ),
        ),
        population=(
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ).with_fitness_value((4, 10, 1)),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ).with_fitness_value((9, 20, 2)),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ).with_fitness_value((16, 30, 4)),
            stk.MoleculeRecord(topology_graph),
        ),
        normalized=(
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ).with_fitness_value(28),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ).with_fitness_value(57),
            stk.MoleculeRecord(
                topology_graph=topology_graph,
            ).with_fitness_value(92),
            stk.MoleculeRecord(topology_graph),
        ),
    )


@pytest.fixture(
    scope='session',
    params=(
        _get_case_data_1,
    ),
)
def sequence(request) -> CaseData:
    return request.param()
