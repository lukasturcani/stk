import numpy as np
import stk


def test_put_caching(mongo_client):
    database_name = "_test_put_caching"
    mongo_client.drop_database(database_name)

    database = stk.ConstructedMoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
    )
    molecule = stk.BuildingBlock("BrCCCBr", [stk.BromoFactory()])
    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=(molecule,),
            repeating_unit="A",
            num_repeating_units=3,
        ),
    )
    database.put(polymer)
    database.put(polymer)

    cache_info = database._put.cache_info()
    assert cache_info.hits == 1
    assert cache_info.misses == 1

    database.put(
        molecule=polymer.with_position_matrix(
            position_matrix=np.zeros((polymer.get_num_atoms(), 3)),
        ),
    )
    cache_info = database._put.cache_info()
    assert cache_info.hits == 1
    assert cache_info.misses == 2


def test_get_caching(mongo_client):
    database_name = "_test_get_caching"
    mongo_client.drop_database(database_name)

    database = stk.ConstructedMoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
    )
    molecule = stk.BuildingBlock("BrCCCBr", [stk.BromoFactory()])
    polymer = stk.ConstructedMolecule(
        topology_graph=stk.polymer.Linear(
            building_blocks=(molecule,),
            repeating_unit="A",
            num_repeating_units=3,
        ),
    )
    database.put(polymer)
    database.get(
        {
            stk.InchiKey().get_key_name(): stk.InchiKey().get_key(polymer),
        }
    )
    database.get(
        {
            stk.InchiKey().get_key_name(): stk.InchiKey().get_key(polymer),
        }
    )

    cache_info = database._get.cache_info()
    assert cache_info.hits == 1
    assert cache_info.misses == 1
