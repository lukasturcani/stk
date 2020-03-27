import pytest
import stk

from ..case_data import CaseData
from ...utilities import MockMongoClient


@pytest.fixture(
    params=(
        CaseData(
            database=stk.ConstructedMoleculeValueMongoDb(
                mongo_client=MockMongoClient(),
                collection='values',
                lru_cache_size=128,
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
        CaseData(
            database=stk.ConstructedMoleculeValueMongoDb(
                mongo_client=MockMongoClient(),
                collection='values',
                lru_cache_size=0,
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
def constructed_molecule_value_mongo_db(request):
    return request.param
