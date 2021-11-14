import pandas as pd
import pytest

import stk

from .case_data import CaseData


def _get_topology_graph() -> stk.polymer.Linear:
    return stk.polymer.Linear(
        building_blocks=(
            stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
        ),
        repeating_unit='A',
        num_repeating_units=2,
    )


def get_generation(*fitness_values):
    v1, v2, v3, *_ = fitness_values
    topology_graph = _get_topology_graph()
    return (
        stk.MoleculeRecord(
            topology_graph=topology_graph,
        ).with_fitness_value(v1),
        stk.MoleculeRecord(
            topology_graph=topology_graph,
        ).with_fitness_value(v2),
        stk.MoleculeRecord(
            topology_graph=topology_graph,
        ).with_fitness_value(v3),
        stk.MoleculeRecord(
            topology_graph=topology_graph,
        )
    )


@pytest.fixture(
    scope='session',
    params=(
        lambda: CaseData(
            plotter=stk.ProgressPlotter(
                generations=(
                    stk.Generation(
                        molecule_records=get_generation(0, 1, 2),
                        mutation_records=(),
                        crossover_records=(),
                    ),
                    stk.Generation(
                        molecule_records=get_generation(10, 20, 30),
                        mutation_records=(),
                        crossover_records=(),
                    ),
                    stk.Generation(
                        molecule_records=get_generation(40, 50, 60),
                        mutation_records=(),
                        crossover_records=(),
                    ),
                    stk.Generation(
                        molecule_records=get_generation(40, 50, 60),
                        mutation_records=(),
                        crossover_records=(),
                    ),
                    stk.Generation(
                        molecule_records=get_generation(70, 80, 90),
                        mutation_records=(),
                        crossover_records=(),
                    ),
                ),
                get_property=lambda record: record.get_fitness_value(),
                y_label='Fitness Value',
                filter=lambda record:
                    record.get_fitness_value() is not None,
            ),
            plot_data=pd.DataFrame({
                'Generation': [0]*3 + [1]*3 + [2]*3 + [3]*3 + [4]*3,
                'Fitness Value': [
                    2., 1., 0.,
                    30., 20., 10.,
                    60., 50., 40.,
                    60., 50., 40.,
                    90., 80., 70.,
                ],
                'Type': ['Max', 'Mean', 'Min']*5
            }),
        ),
    ),
)
def case_data(request) -> CaseData:
    return request.param()
