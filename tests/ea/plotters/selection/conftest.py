import pytest
import stk

from .case_data import CaseData


def get_topology_graph(num_repeating_units):
    return stk.polymer.Linear(
        building_blocks=(stk.BuildingBlock("BrCCBr", [stk.BromoFactory()]),),
        repeating_unit="A",
        num_repeating_units=num_repeating_units,
    )


@pytest.fixture(
    scope="session",
    params=(
        lambda: CaseData(
            selector=stk.Roulette(num_batches=50),
            population=(
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(2),
                ).with_fitness_value(1),
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(3),
                ).with_fitness_value(2),
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(4),
                ).with_fitness_value(3),
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(5),
                ).with_fitness_value(4),
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(6),
                ).with_fitness_value(5),
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(7),
                ).with_fitness_value(6),
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(8),
                ).with_fitness_value(7),
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(9),
                ).with_fitness_value(8),
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(10),
                ).with_fitness_value(9),
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(11),
                ).with_fitness_value(10),
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(12),
                ).with_fitness_value(11),
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(13),
                ).with_fitness_value(12),
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(14),
                ).with_fitness_value(13),
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(15),
                ).with_fitness_value(14),
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(16),
                ).with_fitness_value(15),
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(17),
                ).with_fitness_value(16),
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(18),
                ).with_fitness_value(17),
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(19),
                ).with_fitness_value(18),
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(20),
                ).with_fitness_value(19),
                stk.MoleculeRecord(
                    topology_graph=get_topology_graph(21),
                ).with_fitness_value(20),
            ),
        ),
    ),
)
def case_data(request) -> CaseData:
    return request.param()
