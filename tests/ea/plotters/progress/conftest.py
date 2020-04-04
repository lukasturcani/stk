import pytest
import stk
import pandas as pd

from .case_data import CaseData


@pytest.fixture(
    CaseData(
        plotter=stk.ProgressPlotter(
            generations=(
                stk.Generation(
                    id=0,
                    molecule_records=(
                        stk.MoleculeRecord(
                            molecule=stk.BuildingBlock('BrBr'),
                        ).with_fitness_value(0),
                        stk.MoleculeRecord(
                            molecule=stk.BuildingBlock('BrCBr'),
                        ).with_fitness_value(1),
                        stk.MoleculeRecord(
                            molecule=stk.BuildingBlock('BrCCBr'),
                        ).with_fitness_value(2),
                        stk.MoleculeRecord(
                            molecule=stk.BuildingBlock('BrCCCBr'),
                        ),
                    ),
                ),
                stk.Generation(
                    id=1,
                    molecule_records=(
                        stk.MoleculeRecord(
                            molecule=stk.BuildingBlock('BrBr'),
                        ).with_fitness_value(10),
                        stk.MoleculeRecord(
                            molecule=stk.BuildingBlock('BrCBr'),
                        ).with_fitness_value(20),
                        stk.MoleculeRecord(
                            molecule=stk.BuildingBlock('BrCCBr'),
                        ).with_fitness_value(30),
                        stk.MoleculeRecord(
                            molecule=stk.BuildingBlock('BrCCCBr'),
                        ),
                    ),
                ),
                stk.Generation(
                    id=2,
                    molecule_records=(
                        stk.MoleculeRecord(
                            molecule=stk.BuildingBlock('BrBr'),
                        ).with_fitness_value(10),
                        stk.MoleculeRecord(
                            molecule=stk.BuildingBlock('BrCBr'),
                        ).with_fitness_value(20),
                        stk.MoleculeRecord(
                            molecule=stk.BuildingBlock('BrCCBr'),
                        ).with_fitness_value(30),
                        stk.MoleculeRecord(
                            molecule=stk.BuildingBlock('BrCCCBr'),
                        ),
                    ),
                ),
                stk.Generation(
                    id=3,
                    molecule_records=(
                        stk.MoleculeRecord(
                            molecule=stk.BuildingBlock('BrBr'),
                        ).with_fitness_value(40),
                        stk.MoleculeRecord(
                            molecule=stk.BuildingBlock('BrCBr'),
                        ).with_fitness_value(50),
                        stk.MoleculeRecord(
                            molecule=stk.BuildingBlock('BrCCBr'),
                        ).with_fitness_value(60),
                    ),
                ),
                stk.Generation(
                    id=4,
                    molecule_records=(
                        stk.MoleculeRecord(
                            molecule=stk.BuildingBlock('BrBr'),
                        ).with_fitness_value(70),
                        stk.MoleculeRecord(
                            molecule=stk.BuildingBlock('BrCBr'),
                        ).with_fitness_value(80),
                        stk.MoleculeRecord(
                            molecule=stk.BuildingBlock('BrCCBr'),
                        ).with_fitness_value(90),
                    ),
                )
            ),
            get_property=lambda record: record.get_fitness_value(),
            y_label='Fitness Value',
            filter=lambda record:
                record.get_fitness_value() is not None,
        ),
        plot_data=pd.DataFrame
    )
)
def case_data(request):
    return request.param
