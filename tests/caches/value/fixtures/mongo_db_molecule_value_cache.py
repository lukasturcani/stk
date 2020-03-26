import pytest
import stk

from ..case_data import CaseData
from ...utilities import MockMongoClient


@pytest.fixture(
    params=(
        CaseData(
            cache=stk.MongoDbMoleculeValueCache(
                mongo_client=MockMongoClient(),
                collection='values',
                lru_cache_size=0,
            ),
            molecule=stk.BuildingBlock('BrCCBr'),
            value=12,
        ),
        CaseData(
            cache=stk.MongoDbMoleculeValueCache(
                mongo_client=MockMongoClient(),
                collection='values',
                lru_cache_size=128,
            ),
            molecule=stk.BuildingBlock('BrCCBr'),
            value=12,
        ),
    ),
)
def mongo_db_molecule_value_cache(request):
    return request.param