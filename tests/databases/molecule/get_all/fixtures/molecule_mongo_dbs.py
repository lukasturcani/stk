import pytest
import stk

from ..case_data import CaseData


@pytest.fixture
def molecules():
    return (
        stk.BuildingBlock('CCC'),
        stk.BuildingBlock('BrCCCBr'),
        stk.BuildingBlock('NCCN'),
        stk.BuildingBlock('CCCCC'),
        stk.BuildingBlock('NCCCCN'),
        stk.BuildingBlock('NCC(CCBr)CCN'),
        stk.BuildingBlock('NCCNCCC(CCBr)CCN'),
    )


def get_database(database_name, mongo_client, keys, indices):

    return stk.MoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        jsonizer=stk.MoleculeJsonizer(keys),
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        indices=indices,
    )


@pytest.fixture
def molecule_mongo_dbs(mongo_client, molecules):

    inchi = stk.Inchi()
    smiles = stk.Smiles()
    database_name = '_test_get_entries_molecule'
    mongo_client.drop_database(database_name)

    return CaseData(
        inchi_database=get_database(
            database_name=database_name,
            mongo_client=mongo_client,
            keys=(inchi, ),
            indices=(inchi.get_key_name(), ),
        ),
        smiles_database=get_database(
            database_name=database_name,
            mongo_client=mongo_client,
            keys=(smiles, ),
            indices=(smiles.get_key_name(), ),
        ),
        inchi_key_database=get_database(
            database_name=database_name,
            mongo_client=mongo_client,
            keys=(stk.InchiKey(), ),
            indices=(),
        ),
        inchi_and_smiles_database=get_database(
            database_name=database_name,
            mongo_client=mongo_client,
            keys=(inchi, smiles),
            indices=(),
        ),
        molecules=molecules,
    )
