import numpy as np
import stk


def test_put_caching(mongo_client):
    database_name = "_test_put_caching"
    mongo_client.drop_database(database_name)

    database = stk.MoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        jsonizer=stk.MoleculeJsonizer(
            key_makers=(stk.InchiKey(),),
        ),
    )
    molecule = stk.BuildingBlock("CCC")
    database.put(molecule)
    database.put(molecule)

    cache_info = database._put.cache_info()
    assert cache_info.hits == 1
    assert cache_info.misses == 1

    database.put(
        molecule=molecule.with_position_matrix(
            position_matrix=np.zeros((molecule.get_num_atoms(), 3)),
        ),
    )
    cache_info = database._put.cache_info()
    assert cache_info.hits == 1
    assert cache_info.misses == 2


def test_get_caching(mongo_client):
    database_name = "_test_get_caching"
    mongo_client.drop_database(database_name)

    database = stk.MoleculeMongoDb(
        mongo_client=mongo_client,
        database=database_name,
        jsonizer=stk.MoleculeJsonizer(
            key_makers=(stk.InchiKey(),),
        ),
    )
    molecule = stk.BuildingBlock("CCC")
    database.put(molecule)
    database.get(
        {
            stk.InchiKey().get_key_name(): stk.InchiKey().get_key(molecule),
        }
    )
    database.get(
        {
            stk.InchiKey().get_key_name(): stk.InchiKey().get_key(molecule),
        }
    )

    cache_info = database._get.cache_info()
    assert cache_info.hits == 1
    assert cache_info.misses == 1
