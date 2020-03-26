import pytest
import stk

from ..case_data import CaseData
from ...utilities import MockMongoClient


@pytest.fixture(
    params=(
        CaseData(
            cache=stk.MongoDbConstructedMoleculeValueCache(
                client=MockMongoClient(),
                collection='values',
            ),
            molecule=stk.ConstructedMolecule(
                topology_graph=stk.polymer.Linear(
                    building_blocks=(
                        stk.BuildingBlock(
                            smiles='BrCCBr',
                            functional_groups=[stk.BromoFactory()],
                        ),
                    ),
                    repeating_unit='A',
                    num_repeating_units=2,
                ),
            ),
            value=12,
        ),
    ),
)
def mongo_db_constructed_molecule_value_cache(request):
    return request.param
