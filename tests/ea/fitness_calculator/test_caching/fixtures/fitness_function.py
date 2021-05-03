import pytest
import pymongo
import stk

from .utilities import Counter
from ..case_data import CaseData


_counter = Counter()


def _get_case_data():
    """
    Get a :class:`.CaseData` instance.

    """

    # The basic idea here is that the _counter.get_count method will
    # return a different "fitness value" each time it is called.
    # When the test runs fitness_calculator.get_fitness_value(), if
    # caching is working, the same number as before will be returned.
    # However, if caching is not working, a different number will be
    # returned as the fitness value.

    db = stk.ValueMongoDb(
        mongo_client=pymongo.MongoClient(),
        collection='test_caching',
        database='_stk_pytest_database',
    )

    fitness_calculator = stk.FitnessFunction(
        fitness_function=_counter.get_count,
        input_database=db,
        output_database=db,
    )
    molecule = stk.BuildingBlock('BrCCBr')
    fitness_value = fitness_calculator.get_fitness_value(molecule)

    return CaseData(
        fitness_calculator=fitness_calculator,
        molecule=molecule,
        fitness_value=fitness_value,
    )


@pytest.fixture(
    params=(
        _get_case_data(),
    ),
)
def fitness_function(request):
    return request.param
