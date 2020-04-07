import pytest
import stk

from .case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            selector=stk.Roulette(num_batches=50),
            population=(
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCBr'),
                ).with_fitness_value(10),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCNCBr'),
                ).with_fitness_value(9),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCNNCBr'),
                ).with_fitness_value(2),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCBr'),
                ).with_fitness_value(1),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCCBr'),
                ).with_fitness_value(1),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('ClCBr'),
                ).with_fitness_value(100),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('ClCNCBr'),
                ).with_fitness_value(9),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('ClCNNCBr'),
                ).with_fitness_value(5),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('ClCCBr'),
                ).with_fitness_value(3),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('ClCCCBr'),
                ).with_fitness_value(1),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCl'),
                ).with_fitness_value(50),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCNCCl'),
                ).with_fitness_value(9),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCNNCCl'),
                ).with_fitness_value(2),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCCl'),
                ).with_fitness_value(1),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('BrCCCCl'),
                ).with_fitness_value(1),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('ClCCl'),
                ).with_fitness_value(100),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('ClCNCCl'),
                ).with_fitness_value(9),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('ClCNNCCl'),
                ).with_fitness_value(5),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('ClCCCl'),
                ).with_fitness_value(3),
                stk.MoleculeRecord(
                    molecule=stk.BuildingBlock('ClCCCCl'),
                ).with_fitness_value(1),
            ),
        ),
    ),
)
def case_data(request):
    return request.param
