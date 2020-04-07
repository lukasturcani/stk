import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            selector=stk.Roulette(num_batches=50),
            population=(
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('C'),
                ).with_fitness_value(10),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('C'*2),
                ).with_fitness_value(9),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('C'*3),
                ).with_fitness_value(2),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('C'*4),
                ).with_fitness_value(1),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('C'*5),
                ).with_fitness_value(1),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('C'*6),
                ).with_fitness_value(100),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('C'*7),
                ).with_fitness_value(9),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('C'*8),
                ).with_fitness_value(5),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('C'*9),
                ).with_fitness_value(3),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('C'*10),
                ).with_fitness_value(1),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('C'*11),
                ).with_fitness_value(50),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('C'*12),
                ).with_fitness_value(9),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('C'*13),
                ).with_fitness_value(2),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('C'*14),
                ).with_fitness_value(1),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('C'*15),
                ).with_fitness_value(1),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('C'*16),
                ).with_fitness_value(100),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('C'*17),
                ).with_fitness_value(9),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('C'*18),
                ).with_fitness_value(5),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('C'*19),
                ).with_fitness_value(3),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('C'*20),
                ).with_fitness_value(1),
            ),
        ),
    ),
)
def case_data(request):
    return request.param
