import pytest
import stk

from ..case_data import CaseData


population1 = (
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
)


@pytest.fixture(
    params=(
        CaseData(
            selector=stk.RemoveMolecules(
                remover=stk.Best(2),
                selector=stk.Best(),
            ),
            population=population1,
            selected=(
                stk.Batch(
                    records=(population1[2], ),
                    fitness_values={population1[2]: 2},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[3], ),
                    fitness_values={population1[3]: 1},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[4], ),
                    fitness_values={population1[4]: 1},
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
    ),
)
def remove_molecules(request):
    return request.param
