import pytest
import stk
import rdkit.Chem.AllChem as rdkit

from ...utilities import MockMongoClient
from ..case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            cache=stk.MongoDbMoleculeCache(
                mongo_client=MockMongoClient(),
                lru_cache_size=0,
            ),
            molecule=stk.BuildingBlock('BrCCBr'),
            key={
                'InChIKey':
                    rdkit.MolToInchiKey(rdkit.MolFromSmiles('BrCCBr'))
            },
        ),
        CaseData(
            cache=stk.MongoDbMoleculeCache(
                mongo_client=MockMongoClient(),
                lru_cache_size=0,
                jsonizer=stk.MoleculeJsonizer(
                    key_makers=(
                        stk.MoleculeKeyMaker(
                            key_name='SMILES',
                            get_key=lambda molecule:
                                rdkit.MolToSmiles(
                                    molecule.to_rdkit_mol()
                                )
                        ),
                    ),
                ),
            ),
            molecule=stk.BuildingBlock('BrBr'),
            key={'SMILES': 'BrBr'},
        ),
    ),
)
def mongo_db_molecule_cache(request):
    return request.param
