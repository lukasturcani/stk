from __future__ import annotations

import pytest
import stk

from ..case_data import CaseData
import pymongo


@pytest.fixture
def molecules() -> tuple[stk.ConstructedMolecule]:
    building_blocks = (
        stk.BuildingBlock('BrC#CBr', [stk.BromoFactory()]),
        stk.BuildingBlock('BrCCBr', [stk.BromoFactory()]),
        stk.BuildingBlock('BrCCCBr', [stk.BromoFactory()]),
        stk.BuildingBlock('BrCCCCBr', [stk.BromoFactory()]),
        stk.BuildingBlock('BrCNCBr', [stk.BromoFactory()]),
        stk.BuildingBlock('BrCCNCBr', [stk.BromoFactory()]),
    )

    return [
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(building_blocks[0], ),
                repeating_unit='A',
                num_repeating_units=3,
            ),
        ),
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    building_blocks[0], building_blocks[1]
                ),
                repeating_unit='AB',
                num_repeating_units=3,
            ),
        ),
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(building_blocks[2], ),
                repeating_unit='A',
                num_repeating_units=3,
            ),
        ),
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    building_blocks[2], building_blocks[3]
                ),
                repeating_unit='AB',
                num_repeating_units=3,
            ),
        ),
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(building_blocks[4], ),
                repeating_unit='A',
                num_repeating_units=3,
            ),
        ),
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(
                    building_blocks[4], building_blocks[5]
                ),
                repeating_unit='AB',
                num_repeating_units=3,
            ),
        ),
    ]


def get_database(
    database_name: str,
    mongo_client: pymongo.MongoClient,
    keys: tuple[stk.MoleculeKeyMaker],
    indices: tuple[str],
) -> stk.ConstructedMoleculeMongoDb:

    return stk.ConstructedMoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        jsonizer=stk.ConstructedMoleculeJsonizer(keys),
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        indices=indices,
    )


@pytest.fixture
def constructed_molecule_mongo_dbs(
    mongo_client: pymongo.MongoClient,
    molecules: tuple[stk.ConstructedMolecule],
) -> CaseData:

    inchi = stk.Inchi()
    smiles = stk.Smiles()
    database_name = '_test_get_entries_constructed_molecule'
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
