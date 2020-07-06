import pytest
import stk
import pymongo

from ..case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            database=stk.ValueMongoDb(
                mongo_client=pymongo.MongoClient(),
                collection='values',
                database='_stk_test_database_for_testing',
                put_lru_cache_size=0,
                get_lru_cache_size=0,
            ),
            molecule=stk.BuildingBlock('BrCCBr'),
            value=12,
        ),
        CaseData(
            database=stk.ValueMongoDb(
                mongo_client=pymongo.MongoClient(),
                collection='values',
                database='_stk_test_database_for_testing',
                put_lru_cache_size=128,
                get_lru_cache_size=128,
            ),
            molecule=stk.BuildingBlock('BrCCBr'),
            value=12,
        ),
    ),
)
def mongo_db(request):
    return request.param
