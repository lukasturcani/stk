import pymongo
import pytest
import stk

from ..case_data import CaseData
from .utilities import Counter

_counter = Counter()


def _get_case_data(mongo_client: pymongo.MongoClient) -> CaseData:
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

    fitness_calculator = stk.FitnessFunction(
        fitness_function=_counter.get_count,
        input_database=db,
        output_database=db,
    )
    record = stk.MoleculeRecord(
        topology_graph=stk.polymer.Linear(
            building_blocks=(stk.BuildingBlock("BrCCBr"),),
            repeating_unit="A",
            num_repeating_units=1,
        ),
    )
    fitness_value = fitness_calculator.get_fitness_value(record)

    return CaseData(
        fitness_calculator=fitness_calculator,
        record=record,
        fitness_value=fitness_value,
    )


@pytest.fixture(
    params=(_get_case_data,),
)
def fitness_function(
    request: pytest.FixtureRequest,
    mongo_client: pymongo.MongoClient,
) -> CaseData:
    return request.param(mongo_client)
