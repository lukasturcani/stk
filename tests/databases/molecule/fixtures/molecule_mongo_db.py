import pytest
import stk
import rdkit.Chem.AllChem as rdkit
import pymongo

from ..case_data import CaseData


@pytest.fixture(
    params=(
        CaseData(
            database=stk.MoleculeMongoDb(
                mongo_client=pymongo.MongoClient(),
                database='_test_database',
                put_lru_cache_size=0,
                get_lru_cache_size=0,
            ),
            molecule=stk.BuildingBlock('BrCCBr'),
            key={
                'InChIKey':
                    rdkit.MolToInchiKey(rdkit.MolFromSmiles('BrCCBr'))
            },
        ),
        CaseData(
            database=stk.MoleculeMongoDb(
                mongo_client=pymongo.MongoClient(),
                database='_test_database',
                put_lru_cache_size=0,
                get_lru_cache_size=0,
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
        CaseData(
            database=stk.MoleculeMongoDb(
                mongo_client=pymongo.MongoClient(),
                database='_test_database',
                put_lru_cache_size=128,
                get_lru_cache_size=128,
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
def molecule_mongo_db(request):
    return request.param
