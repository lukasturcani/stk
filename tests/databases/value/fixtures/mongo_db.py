import pytest
import stk

from ..case_data import CaseData
from ...utilities import MockMongoClient


@pytest.fixture(
    params=(
        CaseData(
            database=stk.ValueMongoDb(
                mongo_client=MockMongoClient(),
                collection='values',
                lru_cache_size=0,
            ),
            molecule=stk.BuildingBlock('BrCCBr'),
            value=12,
        ),
        CaseData(
            database=stk.ValueMongoDb(
                mongo_client=MockMongoClient(),
                collection='values',
                lru_cache_size=128,
            ),
            molecule=stk.BuildingBlock('BrCCBr'),
            value=12,
        ),
    ),
)
def molecule_value_mongo_db(request):
    return request.param
