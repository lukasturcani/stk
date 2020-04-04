import pytest
import stk
import pandas as pd

from .case_data import CaseData


def get_generation(*fitness_values):
    v1, v2, v3, *_ = fitness_values
    return (
        stk.MoleculeRecord(
            molecule=stk.BuildingBlock('BrBr'),
        ).with_fitness_value(v1),
        stk.MoleculeRecord(
            molecule=stk.BuildingBlock('BrCBr'),
        ).with_fitness_value(v2),
        stk.MoleculeRecord(
            molecule=stk.BuildingBlock('BrCCBr'),
        ).with_fitness_value(v3),
        stk.MoleculeRecord(
            molecule=stk.BuildingBlock('BrCCCBr'),
        )
    )


@pytest.fixture(
    params=(
        CaseData(
            plotter=stk.ProgressPlotter(
                generations=(
                    stk.Generation(
                        id=0,
                        molecule_records=get_generation(0, 1, 2),
                        mutation_records=(),
                        crossover_records=(),
                    ),
                    stk.Generation(
                        id=1,
                        molecule_records=get_generation(10, 20, 30),
                        mutation_records=(),
                        crossover_records=(),
                    ),
                    stk.Generation(
                        id=2,
                        molecule_records=get_generation(40, 50, 60),
                        mutation_records=(),
                        crossover_records=(),
                    ),
                    stk.Generation(
                        id=3,
                        molecule_records=get_generation(40, 50, 60),
                        mutation_records=(),
                        crossover_records=(),
                    ),
                    stk.Generation(
                        id=4,
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
                'Generation': [0, 1, 2, 3, 4],
                'Fitness Value': [
                    2, 1, 0,
                    30, 20, 10,
                    60, 50, 40,
                    60, 50, 40,
                    90, 80, 70,
                ],
                'Type': ['Max', 'Mean', 'Min']*5
            }),
        )
    ),
)
def case_data(request):
    return request.param
