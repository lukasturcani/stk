import pytest
import stk

from ..case_data import CaseData

topology_graph = stk.polymer.Linear(
    building_blocks=(
        stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
    ),
    repeating_unit='A',
    num_repeating_units=2,
)


@pytest.fixture(
    params=(
        CaseData(
            fitness_normalizer=stk.NormalizerSequence(
                fitness_normalizers=(
                    stk.Multiply((1, 2, 4)),
                    stk.Sum(),
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
                ).with_fitness_value((4, 20, 4)),
                stk.MoleculeRecord(
                    topology_graph=topology_graph,
                ).with_fitness_value((9, 40, 8)),
                stk.MoleculeRecord(
                    topology_graph=topology_graph,
                ).with_fitness_value((16, 60, 16)),
                stk.MoleculeRecord(topology_graph),
            ),
        ),
    ),
)
def sequence(request):
    return request.param
