import pytest
import stk

from ..case_data import CaseData


population1 = (
    stk.MoleculeRecord(
        molecule=stk.BuildingBlock('BrCBr'),
    ).with_fitness_value(5),
    stk.MoleculeRecord(
        molecule=stk.BuildingBlock('BrCNCBr'),
    ).with_fitness_value(3),
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
            selector=stk.AboveAverage(),
            population=population1,
            selected=(
                stk.Batch(
                    records=(population1[0], ),
                    fitness_values={population1[0]: 5},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[0], ),
                    fitness_values={population1[0]: 5},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[1], ),
                    fitness_values={population1[1]: 3},
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        CaseData(
            selector=stk.AboveAverage(num_batches=2),
            population=population1,
            selected=(
                stk.Batch(
                    records=(population1[0], ),
                    fitness_values={population1[0]: 5},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[0], ),
                    fitness_values={population1[0]: 5},
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        CaseData(
            selector=stk.AboveAverage(
                duplicate_molecules=False,
            ),
            population=population1,
            selected=(
                stk.Batch(
                    records=(population1[0], ),
                    fitness_values={population1[0]: 5},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[1], ),
                    fitness_values={population1[1]: 3},
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        CaseData(
            selector=stk.AboveAverage(
                duplicate_batches=False,
            ),
            population=population1,
            selected=(
                stk.Batch(
                    records=(population1[0], ),
                    fitness_values={population1[0]: 5},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[1], ),
                    fitness_values={population1[1]: 3},
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        CaseData(
            selector=stk.AboveAverage(
                batch_size=2,
            ),
            population=population1,
            selected=(
                stk.Batch(
                    records=(population1[0], population1[1]),
                    fitness_values={
                        population1[0]: 5,
                        population1[1]: 3,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[0], population1[2], ),
                    fitness_values={
                        population1[0]: 5,
                        population1[2]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[0], population1[3], ),
                    fitness_values={
                        population1[0]: 5,
                        population1[3]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        CaseData(
            selector=stk.AboveAverage(
                num_batches=2,
                batch_size=2,
            ),
            population=population1,
            selected=(
                stk.Batch(
                    records=(population1[0], population1[1]),
                    fitness_values={
                        population1[0]: 5,
                        population1[1]: 3,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[0], population1[2], ),
                    fitness_values={
                        population1[0]: 5,
                        population1[2]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
    ),
)
def above_average(request):
    return request.param
