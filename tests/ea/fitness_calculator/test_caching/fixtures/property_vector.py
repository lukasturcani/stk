import pytest

import stk

from ..case_data import CaseData
from .utilities import Counter

_counter = Counter()


def _get_case_data(mongo_client):
    """
    Get a :class:`.CaseData` instance.

    Parameters
    ----------
    mongo_client : :class:`pymongo.MongoClient`
        The mongo client the database should connect to.

    """

    # The basic idea here is that the _counter.get_count method will
    # return a different "fitness value" each time it is called.
    # When the test runs fitness_calculator.get_fitness_value(), if
    # caching is working, the same number as before will be returned.
    # However, if caching is not working, a different number will be
    # returned as the fitness value.

    db = stk.ValueMongoDb(
        mongo_client=mongo_client,
        collection="test_caching",
        database="_stk_pytest_database",
    )

    fitness_calculator = stk.PropertyVector(
        property_functions=(_counter.get_count,),
        input_database=db,
        output_database=db,
    )
    molecule = stk.BuildingBlock("BrCCBr")
    fitness_value = fitness_calculator.get_fitness_value(molecule)

    return CaseData(
        fitness_calculator=fitness_calculator,
        molecule=molecule,
        fitness_value=fitness_value,
    )


@pytest.fixture(
    params=(_get_case_data,),
)
def property_vector(request, mongo_client):
    return request.param(mongo_client)
