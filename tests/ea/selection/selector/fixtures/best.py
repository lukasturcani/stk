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
population2 = (
    stk.MoleculeRecord(
        molecule=stk.BuildingBlock('BrCBr'),
    ).with_fitness_value(100),
    stk.MoleculeRecord(
        molecule=stk.BuildingBlock('BrCBr'),
    ).with_fitness_value(100),
    stk.MoleculeRecord(
        molecule=stk.BuildingBlock('BrCNCBr'),
    ).with_fitness_value(1),
)


@pytest.fixture(
    params=(
        CaseData(
            selector=stk.Best(),
            population=population1,
            selected=(
                stk.Batch(
                    records=(population1[0], ),
                    fitness_values={population1[0]: 10},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[1], ),
                    fitness_values={population1[1]: 9},
                    key_maker=stk.Inchi(),
                ),
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
        CaseData(
            selector=stk.Best(num_batches=2),
            population=population1,
            selected=(
                stk.Batch(
                    records=(population1[0], ),
                    fitness_values={population1[0]: 10},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[1], ),
                    fitness_values={population1[1]: 9},
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        CaseData(
            selector=stk.Best(
                duplicate_molecules=False,
            ),
            population=population2,
            selected=(
                stk.Batch(
                    records=(population2[0], ),
                    fitness_values={population2[0]: 100},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population2[2], ),
                    fitness_values={population2[2]: 1},
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        CaseData(
            selector=stk.Best(
                duplicate_batches=False,
            ),
            population=population2,
            selected=(
                stk.Batch(
                    records=(population2[0], ),
                    fitness_values={population2[0]: 100},
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population2[2], ),
                    fitness_values={population2[2]: 1},
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        CaseData(
            selector=stk.Best(
                batch_size=2,
            ),
            population=population1,
            selected=(
                stk.Batch(
                    records=(population1[0], population1[1]),
                    fitness_values={
                        population1[0]: 10,
                        population1[1]: 9,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[0], population1[2], ),
                    fitness_values={
                        population1[0]: 10,
                        population1[2]: 2,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[0], population1[3], ),
                    fitness_values={
                        population1[0]: 10,
                        population1[3]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[0], population1[4], ),
                    fitness_values={
                        population1[0]: 10,
                        population1[4]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[1], population1[2], ),
                    fitness_values={
                        population1[1]: 9,
                        population1[2]: 2,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[1], population1[3], ),
                    fitness_values={
                        population1[1]: 9,
                        population1[3]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[1], population1[4], ),
                    fitness_values={
                        population1[1]: 9,
                        population1[4]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[2], population1[3], ),
                    fitness_values={
                        population1[2]: 2,
                        population1[3]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[2], population1[4], ),
                    fitness_values={
                        population1[2]: 2,
                        population1[4]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[3], population1[4], ),
                    fitness_values={
                        population1[3]: 1,
                        population1[4]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        CaseData(
            selector=stk.Best(
                num_batches=3,
                batch_size=2,
            ),
            population=population1,
            selected=(
                stk.Batch(
                    records=(population1[0], population1[1]),
                    fitness_values={
                        population1[0]: 10,
                        population1[1]: 9,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[0], population1[2], ),
                    fitness_values={
                        population1[0]: 10,
                        population1[2]: 2,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[0], population1[3], ),
                    fitness_values={
                        population1[0]: 10,
                        population1[3]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        CaseData(
            selector=stk.Best(
                batch_size=2,
                duplicate_molecules=False,
            ),
            population=population1,
            selected=(
                stk.Batch(
                    records=(population1[0], population1[1]),
                    fitness_values={
                        population1[0]: 10,
                        population1[1]: 9,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population1[2], population1[3], ),
                    fitness_values={
                        population1[2]: 2,
                        population1[3]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
        CaseData(
            selector=stk.Best(
                batch_size=2,
                duplicate_batches=False,
            ),
            population=population2,
            selected=(
                stk.Batch(
                    records=(population2[0], population2[1]),
                    fitness_values={
                        population2[0]: 100,
                        population2[1]: 100,
                    },
                    key_maker=stk.Inchi(),
                ),
                stk.Batch(
                    records=(population2[0], population2[2], ),
                    fitness_values={
                        population2[0]: 100,
                        population2[2]: 1,
                    },
                    key_maker=stk.Inchi(),
                ),
            ),
        ),
    ),
)
def above_average(request):
    return request.param
