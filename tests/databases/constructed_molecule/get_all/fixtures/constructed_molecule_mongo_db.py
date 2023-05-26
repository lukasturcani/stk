from __future__ import annotations

import pymongo
import pytest
import stk

from ..case_data import CaseData


@pytest.fixture
def molecules() -> tuple[stk.ConstructedMolecule, ...]:
    bb1 = stk.BuildingBlock("BrC#CBr", [stk.BromoFactory()])
    bb2 = stk.BuildingBlock("BrCCBr", [stk.BromoFactory()])
    bb3 = stk.BuildingBlock("BrCCCBr", [stk.BromoFactory()])
    bb4 = stk.BuildingBlock("BrCCCCBr", [stk.BromoFactory()])
    bb5 = stk.BuildingBlock("BrCNCBr", [stk.BromoFactory()])
    bb6 = stk.BuildingBlock("BrCCNCBr", [stk.BromoFactory()])

    return (
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(bb1,),
                repeating_unit="A",
                num_repeating_units=3,
            ),
        ),
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(bb1, bb2),
                repeating_unit="AB",
                num_repeating_units=3,
            ),
        ),
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(bb3,),
                repeating_unit="A",
                num_repeating_units=3,
            ),
        ),
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(bb3, bb4),
                repeating_unit="AB",
                num_repeating_units=3,
            ),
        ),
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(bb5,),
                repeating_unit="A",
                num_repeating_units=3,
            ),
        ),
        stk.ConstructedMolecule(
            topology_graph=stk.polymer.Linear(
                building_blocks=(bb5, bb6),
                repeating_unit="AB",
                num_repeating_units=3,
            ),
        ),
    )


def get_database(
    database_name: str,
    mongo_client: pymongo.MongoClient,
    key_makers: tuple[stk.MoleculeKeyMaker, ...],
    indices: tuple[str, ...],
) -> stk.ConstructedMoleculeMongoDb:
    return stk.ConstructedMoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        jsonizer=stk.ConstructedMoleculeJsonizer(key_makers),
        put_lru_cache_size=0,
        get_lru_cache_size=0,
        indices=indices,
    )


@pytest.fixture
def constructed_molecule_mongo_db(
    mongo_client: pymongo.MongoClient,
    molecules: tuple[stk.ConstructedMolecule, ...],
) -> CaseData:
    inchi = stk.Inchi()
    smiles = stk.Smiles()
    database_name = "_test_get_all_constructed_molecules"
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
