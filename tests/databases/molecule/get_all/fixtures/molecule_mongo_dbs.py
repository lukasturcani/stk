from __future__ import annotations

import pytest
import stk

from ..case_data import CaseData
import pymongo


@pytest.fixture
def molecules() -> tuple[stk.BuildingBlock]:
    return (
        stk.BuildingBlock('CCC'),
        stk.BuildingBlock('BrCCCBr'),
        stk.BuildingBlock('NCCN'),
        stk.BuildingBlock('CCCCC'),
        stk.BuildingBlock('NCCCCN'),
        stk.BuildingBlock('NCC(CCBr)CCN'),
        stk.BuildingBlock('NCCNCCC(CCBr)CCN'),
    )


def get_database(
    database_name: str,
    mongo_client: pymongo.MongoClient,
    keys: tuple[stk.MoleculeKeyMaker],
    indices: tuple[str],
) -> stk.MoleculeMongoDb:

    return stk.MoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        jsonizer=stk.MoleculeJsonizer(keys),
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        indices=indices,
    )


@pytest.fixture
def molecule_mongo_dbs(
    mongo_client: pymongo.MongoClient,
    molecules: tuple[stk.BuildingBlock],
) -> CaseData:

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
