import pytest
import stk
import rdkit.Chem.AllChem as rdkit

from .utilities import MockMongoClient
from ..case_data import CaseData


@pytest.fixture(
    params=(
    ),
)
def mongo_db_constructed_molecule_cache(request):
    return request.param
