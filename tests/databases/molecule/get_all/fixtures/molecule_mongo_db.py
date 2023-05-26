from __future__ import annotations

import pymongo
import pytest
import stk

from ..case_data import CaseData


@pytest.fixture
def molecules() -> tuple[stk.BuildingBlock, ...]:
    return (
        stk.BuildingBlock("CCC"),
        stk.BuildingBlock("BrCCCBr"),
        stk.BuildingBlock("NCCN"),
        stk.BuildingBlock("CCCCC"),
        stk.BuildingBlock("NCCCCN"),
        stk.BuildingBlock("NCC(CCBr)CCN"),
        stk.BuildingBlock("NCCNCCC(CCBr)CCN"),
    )


def get_database(
    database_name: str,
    mongo_client: pymongo.MongoClient,
    key_makers: tuple[stk.MoleculeKeyMaker, ...],
    indices: tuple[str, ...],
) -> stk.MoleculeMongoDb:
    return stk.MoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        jsonizer=stk.MoleculeJsonizer(key_makers),
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        indices=indices,
    )


@pytest.fixture
def molecule_mongo_db(
    mongo_client: pymongo.MongoClient,
    molecules: tuple[stk.BuildingBlock, ...],
) -> CaseData:
    inchi = stk.Inchi()
    smiles = stk.Smiles()
    database_name = "_test_get_all_molecules"
    mongo_client.drop_database(database_name)

    inchi_molecules = molecules[:2]
    smiles_molecules = molecules[2:4]
    inchi_and_smiles_molecules = molecules[4:]

    inchi_database = get_database(
        database_name=database_name,
        mongo_client=mongo_client,
        key_makers=(inchi,),
        indices=(inchi.get_key_name(),),
    )
    smiles_database = get_database(
        database_name=database_name,
        mongo_client=mongo_client,
        key_makers=(smiles,),
        indices=(smiles.get_key_name(),),
    )

    inchi_and_smiles_database = get_database(
        database_name=database_name,
        mongo_client=mongo_client,
        key_makers=(inchi, smiles),
        indices=(),
    )

    for molecule in inchi_molecules:
        inchi_database.put(molecule)

    for molecule in smiles_molecules:
        smiles_database.put(molecule)

    for molecule in inchi_and_smiles_molecules:
        inchi_and_smiles_database.put(molecule)

    inchi_key_database = get_database(
        database_name=database_name,
        mongo_client=mongo_client,
        key_makers=(stk.InchiKey(),),
        indices=(),
    )

    expected_molecules = {
        smiles.get_key(molecule): molecule for molecule in molecules
    }

    return CaseData(
        database=inchi_key_database,
        expected_molecules=expected_molecules,
    )
