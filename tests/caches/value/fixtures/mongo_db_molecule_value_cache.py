import pytest
import stk

from ..case_data import CaseData
from ...utilities import MockMongoClient


@pytest.fixture(
    params=(
        CaseData(
            cache=stk.MongoDbMoleculeValueCache(
                client=MockMongoClient(),
                collection='values',
            ),
            molecule=stk.BuildingBlock('BrCCBr'),
            value=12,
        ),
    ),
)
def mongo_db_molecule_value_cache(request):
    return request.param
